
library(invgamma)
library(fda)
library(MASS)
library(tictoc)
library(pbmcapply)

#' 
#' Bayesian Nonparamteric clustering of functional data.
#'
#' @param n_iter number of iterations
#' @param burnin number of burnin iterations 
#' @param thin thinning step
#' @param M stick-breaking representation
#' @param mass total mass parameter
#' 
#' @param smoothing a list with the following components:
#'                  X data matrix n x T
#'                  basis basisfd object from fda library 
#'                  beta coefficients of projection in basis of X
#'                  time.grid vector defining the time grid
#' 
#' @param hyperparam a list with hyperparameters specification
#' 
#' 
#' @return a list with the following components:
#'         K: (n_iter-burnin) x n matrix, [K]ij = cluster of observation j at iteration i
#'         algorithm_parameters
#' 

FBNP_hyper <- function (n_iter, 
                  burnin=0, 
                  thin=1, 
                  M, 
                  mass,
                  smoothing,
                  hyperparam)
  
{
  
  #### DATA ---------------------------------------------------------------------------------------
  
  X <- smoothing$X
  basis <- smoothing$basis
  beta <- smoothing$beta
  time.grid <- smoothing$time.grid
  
  ##### PARAMETERS SETTING ------------------------------------------------------------------------
  
  # number of data
  n <- dim(X)[1]
  # time steps
  n_time <- dim(X)[2]
  # number of basis functions
  L <- basis$nbasis
  
  ## evaluate the basis on the time grid
  # basis.t: (L x n_time) matrix of basis fuctions evaluated on the time time grid
  basis.t <- t(eval.basis(time.grid, basis))
  
  #### HYPERPARAMETERS ----------------------------------------------------------------------------
  
  a           <- hyperparam$a
  b           <- hyperparam$b
  c           <- hyperparam$c 
  d           <- hyperparam$d
  #m0          <- hyperparam$m0
  #Lambda0     <- hyperparam$Lambda0
  #Lambda0_inv <- solve(Lambda0)
  theta0 <- hyperparam$m0
  k0     <- L
  nu0    <- L
  delta0 <- hyperparam$Lambda0
  
  #### LATENT RV'S INITIALIZATION -----------------------------------------------------------------
  
  ## SIGMA2
  # numeric(M)
  # (sigma2)_i: sigma^2 del cluster i-esimo
  sigma2 <- rinvgamma(n = M, shape=a, rate=b)
  
  ## PHI
  # M x n_time matrix
  # (phi)_ij: phi_t del cluster i-esimo al time point j-esimo (j=1:1600)
  phi <- matrix(rinvgamma(n=M*n_time, shape=c, rate=d), byrow=TRUE, nrow=M, ncol=n_time)
  
  ## mu_coef (coefficients of basis projection of mu)
  mu_coef <- matrix(0, nrow=M, ncol=L) # sample coefficients of basis projection
  for(j in 1:M)
  {
    # sample Lambda0 (covariance of mu_coef[j,])
    Lambda0 <- LaplacesDemon::rinvwishart(nu0, delta0)
    # sample m0 (mean of mu_coef[j,])
    m0 <- mvrnorm(n=1, mu=theta0, Sigma=Lambda0)
    # sample mu_coef[j,] (coefficients of basis projection of mu)
    mu_coef[j,] <- mvrnorm(n=1, mu=m0, Sigma=Lambda0)
  }
  
  #plot(time.grid, mu_coef[1,]%*%basis.t, type='l')
  #plot(time.grid, t(basis.t)%*%mu_coef[1,], type='l')
  
  
  ## MU
  # (mu)it = mu_i(t)
  # M x L matrix
  # evaluation of mu on the time grid
  
  # sbagliato:
  #mu <- matrix(mu_coef %*% basis.t, byrow=TRUE, nrow=M, ncol=n_time) 
  
  # giusto:
  mu <- matrix(mu_coef %*% basis.t, byrow=FALSE, nrow=M, ncol=n_time)
  #plot(time.grid, mu[1,], type='l')
  
  #### ALGORITHM PARAMETERS -----------------------------------------------------------------------
  
  # Stick breaking
  V <- numeric(M)
  p <- rep(1/M, M)
  
  #(K)_ij: cluster di assegnazione all'iterazione i dell'osservazione j
  K <- matrix(0, nrow=(n_iter-burnin), ncol=n)
  
  # K_curr salva l'assegnazione corrente, così che possiamo salvare i K solo
  # per le iterazioni dopo il burnin
  K_curr <- sample(1:n)
  #K_curr <- rep(1,n)
  
  #### ALGORITHM ----------------------------------------------------------------------------------
  
  pb <- progressBar(0, max = n_iter, initial = 0, style = "ETA")
  tic(quiet=TRUE)
  
  
  for(iter in 1:n_iter)
  {
    
    #### STEP 1: acceleration step ####------------------------------------------------------------
    
    # iterate over kernels
    for(j in 1:M)
    {
      # observations in cluster j
      indexes_j <- which(K_curr==j)
      # number of observations in cluster j
      r <- length(indexes_j)
      
      # if the kernel is not empty
      if(r>0)
      {
        ## SIGMA2
        a_r <- a + r*n_time*0.5
        b_r <- b
        for(g in indexes_j)
        {
          b_r <- b_r + 0.5*sum( (X[g,] - mu[j,])^2/phi[j,] ) # somma sui tempi
        }
        sigma2[j] <- rinvgamma(n = 1, shape=a_r, rate=b_r)
        
        ## PHI
        c_r <- c + r*0.5
        for(t in 1:n_time)
        {
          d_r <- d + sum( (X[indexes_j,t] - mu[j,t])^2 )/(2*sigma2[j]) # somma su indexes
          phi[j,t] <- rinvgamma(n=1, shape=c_r, rate=d_r)
        }
        
        ## Lambda0
        nu.n <- nu0 + 1
        k.n <- k0 + 1
        theta.n <- (k0*theta0 + mu_coef[j,])/k.n
        delta.n <- delta0 + k0/k.n * (mu_coef[j,]-theta0) %*% t(mu_coef[j,]-theta0)
        Lambda0 <- LaplacesDemon::rinvwishart(nu.n, delta.n) # TODO: change name
        Lambda0_inv <- solve(Lambda0)
        
        ## m0
        m0 <- mvrnorm(n=1, mu=theta.n, Sigma=Lambda0) # TODO: change name
        
        ## mu_coef (coefficients of basis projection)
        # first compute coefficients of the posterior distribution of mu_coef:
        magic <- matrix(0, nrow=L, ncol=L)
        for(t in 1:n_time)
        {
          magic <- magic + ( basis.t[,t] %*% t(basis.t[,t]) )/(phi[j,t])
        }
        magic <- magic/sigma2[j]
        Lambda_r <- solve(Lambda0_inv + r*magic)
        if(r>1)
        {
          mu_r <- Lambda_r %*% ( Lambda0_inv %*% m0 + magic %*% colSums(beta[indexes_j,]) )
        } else {
          mu_r <- Lambda_r %*% ( Lambda0_inv %*% m0 + magic %*% beta[indexes_j,] )
        }
        # then sample mu_coef:
        mu_coef[j,] <- mvrnorm(n=1, mu=mu_r, Sigma=Lambda_r)
        
        ## MU
        mu[j,] <- mu_coef[j,] %*% basis.t # evaluation of mu on the time grid
        
      } else {
        
        ## SIGMA2
        sigma2[j] <- rinvgamma(n = 1, shape=a, rate=b)
        
        ## PHI
        phi[j,] <- rinvgamma(n=n_time, shape=c, rate=d)
        
        ## Lambda0
        Lambda0 <- LaplacesDemon::rinvwishart(nu0, delta0) # TODO:LAMBDA
        
        ## m0
        m0 <- mvrnorm(n=1, mu=theta0, Sigma=Lambda0)
        
        ## mu_coef: sample coefficients of basis projection
        mu_coef[j,] <- mvrnorm(n=1, mu=m0, Sigma=Lambda0)
        
        ## MU
        mu[j,] <- mu_coef[j,] %*% basis.t # evaluation of mu on the time grid
      }
      
    }
    
    #### STEP 2: K ####----------------------------------------------------------------------------
    
    # (pi)j: probability that observation i belongs to cluster j, Ki~discrete(p_i[1],...,p_i[M]) 
    #        (lo sovrascriviamo per ogni i dentro il ciclo)
    p_i <- numeric(M)
    
    # iterate over the observations
    for(i in 1:n)
    {
      # iterate over kernels
      for(j in 1:M)
      {
        p_i[j] <- log(p[j]) + sum( (-0.5)*log(2*pi*sigma2[j]*phi[j,]) - ((X[i,]-mu[j,])^2)/(2*sigma2[j]*phi[j,]) ) # sum over times
      }
      #p_i <- p_i - max(p_i) # TODO: avoid the subtraction of max and the application of exponential
      #p_i <- p_i - min(p_i)
      p_i <- exp(p_i)#/sum(exp(p_i))
      
      # update the cluster assignment at current iteration
      K_curr[i] <- sample(1:M, size=1, prob=p_i)
      
      # save K
      if(iter > burnin)
      {
        K[iter-burnin,i] <- K_curr[i]
      }
      
    }
    
    #### STEP 3: weights ####----------------------------------------------------------------------
    
    V[1] <- rbeta(1, 1 + sum(K_curr == 1), mass + sum(K_curr > 1))
    p[1] <- V[1]
    
    for(l in 2:(M - 1))
    {
      V[l] <- rbeta(1, 1 + sum(K_curr == l), mass + sum(K_curr > l))
      p[l] <- V[l] * prod(1 - V[1:(l - 1)])
    }
    
    V[M] <- 1
    p[M] <- V[M] * prod(1 - V[1:(M - 1)])
    
    # ProgressBar
    setTxtProgressBar(pb, iter)
    
  }
  
  elapsed <- toc(quiet=TRUE)
  print(paste0("Elapsed: ", round((elapsed$toc - elapsed$tic)/60), 
               " min, ", round((elapsed$toc - elapsed$tic)%%60), " sec."))
  
  #### RETURN -------------------------------------------------------------------------------------
  
  algo_parameters <- c('n_iter' = n_iter,
                       'burnin' = burnin,
                       'thinning' = thin,
                       'M' = M,
                       'mass' = mass)
  
  out <- list(K, algo_parameters)
  names(out) <- c("K", "algorithm_parameters")
  
  return(out)
  
}




#' TODO list:
#'             - consider removing thinning parameter
#'             - consider the more general case in which basis is replace with fdParobj,
#'               a generic functional parameter object, in order to allow also roughness penalization smoothing
#'               (in tal caso non possiamo semplicemente valutare i coefficienti in basi moltiplicando matrici)
#' 

# Adesso stiamo salvando in mu_coef_out soltatno mu_coef dell'ultimo cluster valutato
# Quindi devo salvare mu_coef come matrice, che in ogni riga salva i coefficienti di mu del cluster j_esimo
# In questo modo, avrò una lista di matrici

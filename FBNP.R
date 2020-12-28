
library(invgamma)
library(fda)
library(MASS)
library(tictoc)
library(pbmcapply)

#' 
#' Bayesian Nonparamteric clustering of functional data.
#'
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
#'         K:            (n_iter-burnin) x n matrix, [K]ij = cluster of observation j at iteration i
#'         mu_coef_out:  (n_iter-burnin)-long list, keeps track of mu_coef
#'         sigma2_out:   (n_iter-burnin)-long list, keeps track og sigma2
#'         probs_j_out:  (n_iter-burnin)-long list, keeps track of the probability of each cluster
#'         probs_ij_out: (n_iter-burnin)-long list, keeps track of the probabilities (matrix) of each observation to fall in each cluster  
#'         algorithm_parameters
#' 

FBNP <- function (n_iter, burnin=0, thin=1, M, mass,
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
  # basis2.t: n_time vector, b2[t] is the dot product between basis.t[,t] and itself
  basis2.t <- sapply(1:n_time, function(t) crossprod(basis.t[,t], basis.t[,t]))
  
  #### HYPERPARAMETERS ----------------------------------------------------------------------------
  
  a           <- hyperparam$a
  b           <- hyperparam$b
  c           <- hyperparam$c 
  d           <- hyperparam$d
  m0          <- hyperparam$m0
  Lambda0     <- hyperparam$Lambda0
  Lambda0_inv <- solve(Lambda0)
  
  
  #### LATENT RV'S INITIALIZATION -----------------------------------------------------------------
  
  ## SIGMA2
  # numeric(M)
  # (sigma2)_i: sigma^2 del cluster i-esimo
  sigma2 <- rinvgamma(n = M, shape=a, scale=b)
  
  ## PHI
  # M x n_time matrix
  # (phi)_ij: phi_t del cluster i-esimo al time point j-esimo (j=1:1600)
  phi <- matrix(rinvgamma(n=M*n_time, shape=c, scale=d), byrow=TRUE, nrow=M, ncol=n_time)
  
  ## MU
  # M x L matrix
  # (mu)it = mu_i(t)
  mu_coef <- matrix(mvrnorm(n=M, mu=m0, Sigma=Lambda0), byrow=TRUE, nrow=M, ncol=L) # sample coefficients of basis projection
  mu <- matrix(mu_coef %*% basis.t, byrow=TRUE, nrow=M, ncol=n_time) # evaluation of mu on the time grid
  
  #### ALGORITHM PARAMETERS -----------------------------------------------------------------------
  
  # Stick breaking
  V <- numeric(M)
  p <- rep(1/M, M)
  
  #(K)_ij: cluster di assegnazione all'iterazione i dell'osservazione j
  K <- matrix(0,nrow=(n_iter-burnin),ncol=n)
  
  # K_curr salva l'assegnazione corrente, così che possiamo salvare i K solo
  # per le iterazioni dopo il burnin
  K_curr <- sample(1:n) # perché non facciamo sample 1:M?
  
  #### RETURN VARIABLES ---------------------------------------------------------------------------
  
  mu_coef_out <- list()
  sigma2_out  <- list()
  
  # matrix that tells for each observation the probabilities of observation i belonging to cluster j after STEP 2
  probs_ij <- matrix(0, nrow = n, ncol = M)
  
  # a list containing the previous matrix for each iteration
  probs_ij_out <- list()
  
  # a list containing the probabilities p of belonging to a given cluster
  probs_j_out  <- list()
  
  
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
          b_r <- b_r + 0.5*sum( (X[g,] - mu[j,])^2/phi[j,] )
        }
        sigma2[j] <- rinvgamma(n = 1, shape=a_r, scale=b_r)
        
        ## PHI
        c_r <- c + r*n_time*0.5
        for(t in 1:n_time)
        {
          d_r <- d + sum( (X[indexes_j,t] - mu[j,t])^2/(2*sigma2[j]) )
          phi[j,t] <- rinvgamma(n=1, shape=c_r, scale=d_r)
        }
        
        ## MU
        magic <- 1/sigma2[j] * sum( basis2.t/phi[j,] )
        
        Lambda_r <- solve(Lambda0_inv + r*magic*diag(L))
        if(r>1)
        {
          mu_r <- Lambda_r %*% ( Lambda0_inv %*% m0 + magic * colSums(beta[indexes_j,]) )
        } else {
          mu_r <- Lambda_r %*% ( Lambda0_inv %*% m0 + magic * beta[indexes_j,] )
        }
        
        # sample coefficients of basis projection
        mu_coef[j,] <- mvrnorm(n=1, mu=mu_r, Sigma=Lambda_r)
        
        # evaluation of mu on the time grid
        mu[j,] <- mu_coef[j,] %*% basis.t
        
      } else {
        ## SIGMA2
        sigma2[j] <- rinvgamma(n = 1, shape=a, scale=b)
        
        ## PHI
        phi[j,] <- rinvgamma(n=n_time, shape=c, scale=d)
        
        ## MU
        # sample coefficients of basis projection
        mu_coef[j,] <- mvrnorm(n=1, mu=m0, Sigma=Lambda0)
        # evaluation of mu on the time grid
        mu[j,] <- mu_coef[j,] %*% basis.t
      }
        
    }
    
    #### STEP 2: K ####----------------------------------------------------------------------------
    
    # (pi)j: probability that observation i belongs to cluster j
    p_i <- numeric(M)
    
    # iterate over the observations
    for(i in 1:n)
    {
      # iterate over kernels
      for(j in 1:M)
      {
        p_i[j] <- log(p[j]) + sum( (-0.5)*log(2*pi*sigma2[j]*phi[j,]) - (X[i,]-mu[j,])^2/(2*sigma2[j]*phi[j,]) )
      }
      # TODO: avoid the subtraction of max and the application of exponential
      p_i <- p_i - max(p_i) # TODO: controllare
      p_i <- exp(p_i)/sum(exp(p_i))
      # cluster assignment at current iteration
      K_curr[i] <- sample(1:M, size=1, prob=p_i)
      if(iter > burnin)
      {
        K[iter - burnin,i] <- K_curr[i]
      }
      
      
      # save
      if(iter > burnin)
      {
        probs_ij[i,] <- p_i 
      }
      
    }
    
    # save
    if(iter > burnin)
    {
      probs_ij_out[[iter - burnin]] <- probs_ij
    }
    
    #### STEP 3: weights ####----------------------------------------------------------------------
    
    V[1] <- rbeta(1, 1 + sum(K_curr == 1), mass + sum(K_curr > 1))
    p[1] <- V[1]
    
    for(l in 2:(M - 1))
    {
      V[l] <- rbeta(1, 1 + sum(K_curr == l), mass + sum(K_curr > j))
      p[l] <- V[l] * prod(1 - V[1:(l - 1)])
    }
    
    V[M] <- 1
    p[M] <- V[M] * prod(1 - V[1:(M - 1)])
    
    
    # save
    if(iter > burnin)
    {
      mu_coef_out[[iter - burnin]] <- mu_coef
      sigma2_out [[iter - burnin]] <- sigma2
      probs_j_out[[iter - burnin]] <- p
    }
    
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
  
  out <- list(K, mu_coef_out, sigma2_out, probs_j_out, probs_ij_out, algo_parameters)
  names(out) <- c("K", "mu_coef_out", "sigma2_out", "probs_j_out", "probs_ij_out", "algorithm_parameters")
  

  
  return(out)

}




#' TODO list:
#'             - consider removing thinning parameter
#'             - consider the more general case in which basis is replace with fdParobj,
#'               a generic functional parameter object, in order to allow also roughness penalization smoothing
#'               (in tal caso non possiamo semplicemente valutare i coefficienti in basi moltiplicando matrici)
#' 


## PREVIOUS PRIOR ELICITATION:
## hyperparameters of mu
#m0 <- rep(0,L)
#xi <- 1
#Lambda0 <- diag(xi,L)
#Lambda0_inv <- diag(1/xi,L) 
#
## hyperparameters of sigma^2
#a <- 2.01
#b <- 1.01
#
## hyperparameters of phi_t
#c <- 2.01
#d <- 1.01

# Adesso stiamo salvando in mu_coef_out soltatno mu_coef dell'ultimo cluster valutato
# Quindi devo salvare mu_coef come matrice, che in ogni riga salva i coefficienti di mu del cluster j_esimo
# In questo modo, avrò una lista di matrici

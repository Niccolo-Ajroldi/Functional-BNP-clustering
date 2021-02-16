
library(invgamma)
library(fda)
library(MASS)
library(pbmcapply)

#' 
#' Bayesian Nonparamteric clustering of functional data.
#'
#' @param n_iter number of iterations
#' @param burnin number of burnin iterations 
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
#'         K:                    (n_iter-burnin) x n matrix, [K]ij = cluster of observation j at iteration i
#'         log_L:                (n_iter-burnin)-long vector, keeps track of the overall log-likelihood
#'         COUNTER:              (n_iter-burnin)-long vector, keeps track of the number of new proposed clusters at each iteration
#'         algorithm_parameters: list of the input parameters
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
  
  # smooth data:
  X <- beta %*% basis.t
  
  #### HYPERPARAMETERS ----------------------------------------------------------------------------
  
  c           <- hyperparam$c 
  d           <- hyperparam$d
  m0          <- hyperparam$m0
  Lambda0     <- hyperparam$Lambda0
  Lambda0_inv <- solve(Lambda0)
  
  #### LATENT RV'S INITIALIZATION -----------------------------------------------------------------
  
  ## PHI
  # M x n_time matrix
  # (phi)_ij: phi_t del cluster i-esimo al time point j-esimo (j=1:1600)
  phi <- matrix(rinvgamma(n=M*n_time, shape=c, rate=d), byrow=TRUE, nrow=M, ncol=n_time)
  
  ## MU
  # M x L matrix
  # (mu)it = mu_i(t)
  mu_coef <- matrix(mvrnorm(n=M, mu=m0, Sigma=Lambda0), byrow=TRUE, nrow=M, ncol=L) # sample coefficients of basis projection
  mu <- matrix(mu_coef %*% basis.t, byrow=TRUE, nrow=M, ncol=n_time) # evaluation of mu on the time grid
  
  #### ALGORITHM PARAMETERS -----------------------------------------------------------------------
  
  # Stick breaking
  V <- numeric(M)
  p <- rep(1/M, M)
  
  #(K)_ij: cluster assignment at iteration i of observation j
  K <- matrix(0,nrow=(n_iter-burnin),ncol=n)
  
  # save current cluster assignment
  K_curr <- sample(1:M, size=n)
  
  # save next cluster assignment, it is an information used for the evaluation
  # of COUNTER (see RETURN VARIABLES)
  K_new <- K_curr
  
  #### RETURN VARIABLES ---------------------------------------------------------------------------

  # a vector containing the evaluation of the log-likelihood at each iteration
  overall_logL <- numeric(n_iter-burnin)
  
  # a vector containing the number of new clusters (i.e. not seen in the iteration before)
  # proposed at each iteration
  COUNTER <-  numeric(n_iter-burnin)
  
  #### ALGORITHM ----------------------------------------------------------------------------------
  
  pb <- progressBar(0, max = n_iter, initial = 0, style = "ETA")
  
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
        
        ## PHI
        c_r <- c + r*0.5
        for(t in 1:n_time)
        {
          d_r <- d + sum( (X[indexes_j,t] - mu[j,t])^2/2 ) # somma su indexes
          phi[j,t] <- rinvgamma(n=1, shape=c_r, rate=d_r)
        }
        
        ## MU
        magic <- matrix(0, nrow=L, ncol=L)
        for(t in 1:n_time)
        {
          magic <- magic + ( basis.t[,t] %*% t(basis.t[,t]) )/(phi[j,t])
        }
        
        Lambda_r <- solve(Lambda0_inv + r*magic)
        
        if(r>1)
        {
          mu_r <- Lambda_r %*% ( Lambda0_inv %*% m0 + magic %*% colSums(beta[indexes_j,]) )
        } else {
          mu_r <- Lambda_r %*% ( Lambda0_inv %*% m0 + magic %*% beta[indexes_j,] )
        }
        
        # sample coefficients of basis projection
        mu_coef[j,] <- mvrnorm(n=1, mu=mu_r, Sigma=Lambda_r)
        
        # evaluation of mu on the time grid
        mu[j,] <- mu_coef[j,] %*% basis.t
        
      } else {
        
        ## PHI
        phi[j,] <- rinvgamma(n=n_time, shape=c, rate=d)
        
        ## MU
        # sample coefficients of basis projection
        mu_coef[j,] <- mvrnorm(n=1, mu=m0, Sigma=Lambda0)
        # evaluation of mu on the time grid
        mu[j,] <- mu_coef[j,] %*% basis.t
      }
      
    }
    
    #### STEP 2: K ####----------------------------------------------------------------------------
    
    # (pi)j: probability that observation i belongs to cluster j, Ki~discrete(p_i[1],...,p_i[M]) 
    #        (we overwrite it for each observation i in the loop)
    p_i <- numeric(M)
    logL <- 0
    counter <- 0
    
    # iterate over the observations
    for(i in 1:n)
    {
      # iterate over kernels
      for(j in 1:M)
      {
        p_i[j] <- log(p[j]) + sum( (-0.5)*log(2*pi*phi[j,]) - ((X[i,]-mu[j,])^2)/(2*phi[j,]) )
      }
      
      p_i <- exp(p_i)
      
      K_new[i] <- sample(1:M, size=1, prob=p_i)
      
      # update the overall log-likelihood for diagnostic
      logL <- logL + sum( (-0.5)*log(2*pi*phi[K_new[i],]) - ((X[i,]-mu[K_new[i],])^2)/(2*phi[K_new[i],]) )
      
      # save
      if(iter > burnin)
      {
        K[iter-burnin,i] <- K_new[i]
      }
     
    }

    counter <- length(unique(K_new)) - sum(unique(K_new) %in% unique(K_curr))
    K_curr <- K_new
    
    # save
    if(iter > burnin)
    {
      overall_logL[iter-burnin] <- logL
      COUNTER[iter-burnin] <- counter
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

  #### RETURN -------------------------------------------------------------------------------------
  
  algorithm_parameters <- c('n_iter' = n_iter,
                       'burnin' = burnin,
                       'M' = M,
                       'mass' = mass)
  
  out <- list(K, overall_logL, COUNTER, algorithm_parameters)
  names(out) <- c("K", "logL", "counter", "algorithm_parameters")
  
  return(out)
  
}

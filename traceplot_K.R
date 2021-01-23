
traceplot_K <- function (out,
                         smoothing_list,
                         run_parameters
                         )  
{
  
  # get variables
  K            <- out$K
  mu_coef_out  <- out$mu_coef_out
  sigma2_out   <- out$sigma2_out
  probs_j_out  <- out$probs_j_out
  probs_ij_out <- out$probs_ij_out
  phi_out      <- out$phi_out
  
  n         <- dim(smoothing_list$X)[1]
  basis     <- smoothing_list$basis
  time.grid <- smoothing_list$time.grid
  n_time    <- length(time.grid)
  
  n_iter <- run_parameters$algorithm_parameters[[1]]
  burnin <- run_parameters$algorithm_parameters[[2]]
  thin   <- run_parameters$algorithm_parameters[[3]]
  M      <- run_parameters$algorithm_parameters[[4]]
  mass   <- run_parameters$algorithm_parameters[[5]]
  
  # faccio due plot del cluster assignment, così non mi da prolemi di margin too large
  nhalf <- round(n/2)
  
  # first half of observation
  X11(width=1300, height=700)
  par(mfrow=n2mfrow(n/2))
  par(oma=c(0,0,2,0))
  for(i in 1:nhalf)
  {
    traceplot(as.mcmc(K[,i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  }
  title("Cluster assignment", outer = TRUE)
  
  # second half of observations
  X11(width=1300, height=700)
  par(mfrow=n2mfrow(n-nhalf))
  par(oma=c(0,0,2,0))
  for(i in (nhalf+1):n)
  {
    traceplot(as.mcmc(K[,i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  }
  title("Cluster assignment", outer = TRUE)
  
  
}



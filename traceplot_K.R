
traceplot_K <- function (out,
                         smoothing_list,
                         run_parameters,
                         plot.title=""
                         )  
{
  
  # extract content of out
  K         <- out$K
  n         <- dim(smoothing_list$X)[1]
  basis     <- smoothing_list$basis
  time.grid <- smoothing_list$time.grid
  n_time    <- length(time.grid)
  n_iter    <- run_parameters$algorithm_parameters[[1]]
  burnin    <- run_parameters$algorithm_parameters[[2]]
  M         <- run_parameters$algorithm_parameters[[3]]
  mass      <- run_parameters$algorithm_parameters[[4]]
  
  # faccio due plot del cluster assignment, così non mi da prolemi di margin too large
  nhalf <- round(n/2)
  
  # first half of observation
  X11(width=1300, height=700)
  par(mfrow=n2mfrow(n/2))
  par(oma=c(0,0,2,0))
  for(i in 1:nhalf)
    traceplot(as.mcmc(K[-1,i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  #title("Cluster allocation variables", outer = TRUE)
  title(paste0("Cluster allocation variables ", plot.title), outer = TRUE)
  
  # second half of observations
  X11(width=1300, height=700)
  par(mfrow=n2mfrow(n-nhalf))
  par(oma=c(0,0,2,0))
  for(i in (nhalf+1):n)
    traceplot(as.mcmc(K[-1,i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  #title("Cluster allocation variables", outer = TRUE)
  title(paste0("Cluster allocation variables ", plot.title), outer = TRUE)
  
  
}



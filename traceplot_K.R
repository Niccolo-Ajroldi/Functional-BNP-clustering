
#'
#' Traceplot_K
#' 
#' Given the output of the run and the list of smoothing parameters used for the
#' construction of the smoothed data, the function produces traceplots of the cluster
#' assignment of each observation, of the overall log-likelihood and of the COUNTER
#' 
#' @param out: output of the run, contains K, logL, COUNTER, algorithm_parameters
#' @param smoothing_list: contains the output of function smoothing, here we use basis
#'                        and time.grid
#'



traceplot_K <- function (out,
                         smoothing_list,
                         plot.title="")  
{
  
  # extract content of out
  K         <- out$K
  logL      <- out$logL
  counter   <- out$counter
  n         <- dim(smoothing_list$X)[1]
  basis     <- smoothing_list$basis
  time.grid <- smoothing_list$time.grid
  n_time    <- length(time.grid)
  n_iter    <- out$algorithm_parameters[[1]]
  burnin    <- out$algorithm_parameters[[2]]
  M         <- out$algorithm_parameters[[3]]
  mass      <- out$algorithm_parameters[[4]]
  
  # faccio due plot del cluster assignment, così non mi da prolemi di margin too large
  nhalf <- round(n/2)
  
  # first sixth of observation
  X11(width=1300, height=700)
  par(mfrow=n2mfrow(nhalf))
  par(oma=c(0,0,2,0))
  for(i in 1:nhalf)
    traceplot(as.mcmc(K[-20,i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  #title("Cluster allocation variables", outer = TRUE)
  title(paste0("Cluster allocation variables ", plot.title), outer = TRUE)
  
  # second half of observations
  X11(width=1300, height=700)
  par(mfrow=n2mfrow(nhalf))
  par(oma=c(0,0,2,0))
  for(i in (nhalf+1):(2*nhalf))
    traceplot(as.mcmc(K[-1,i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  #title("Cluster allocation variables", outer = TRUE)
  title(paste0("Cluster allocation variables ", plot.title), outer = TRUE)
  
  x11()
  par(mfrow=c(2,1))
  traceplot(as.mcmc(logL), main="Traceplot for the logLikelihood")
  traceplot(as.mcmc(counter), main="Traceplot of counter")
  
}



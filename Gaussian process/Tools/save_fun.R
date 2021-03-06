#'
#' Save Function
#'
#' It is a utility function that helps in saving the results of the
#' run, along with the input parameters, the PSM and some traceplots useful 
#' for a quick diagnostic
#'
#' @param out: output of the run of the model
#' @param name_dir: directory in which to save
#' 


save_fun <- function(out, name_dir)
{
  source("Posterior inference/Tools/PSM.R")
  library(coda)
  library(devtools)
  library(mcclust.ext)
  
  # save old directory to get back there at the end
  dir.current <- getwd()
  
  # name of directory where I will put plots, I use current time in the name
  new.dir <- paste0(dir.current,name_dir)
  
  # create such directory and go there
  dir.create(new.dir)
  setwd(new.dir) 
  
  # save RData ----------------------------------------------------------------------
  algorithm_parameters <- out$algorithm_parameters
  save(out, algorithm_parameters, file = "Output.RData")
  
  # save traceplots -----------------------------------------------------------------
  # extract content of out
  K         <- out$K
  n         <- dim(smoothing_list$X)[1]
  basis     <- smoothing_list$basis
  time.grid <- smoothing_list$time.grid
  n_time    <- length(time.grid)
  n_iter    <- out$algorithm_parameters[[1]]
  burnin    <- out$algorithm_parameters[[2]]
  M         <- out$algorithm_parameters[[3]]
  mass      <- out$algorithm_parameters[[4]]
  
  
  # faccio due plot del cluster assignment, cos� non mi da prolemi di margin too large
  nhalf <- round(n/2)
  
  # first half of observation
  png(file = "Traceplot_K_1.png", width = 8000, height = 5000, units = "px", res = 800)
  par(mfrow=n2mfrow(n/2))
  par(oma=c(0,0,2,0))
  for(i in 1:nhalf)
    traceplot(as.mcmc(K[,i]), main=paste0("Observation ",i))

  title(paste0("Cluster allocation variables "), outer = TRUE)
  dev.off()
  
  # second half of observations
  png(file = "Traceplot_K_2.png", width = 8000, height = 5000, units = "px", res = 800)
  par(mfrow=n2mfrow(n-nhalf))
  par(oma=c(0,0,2,0))
  for(i in (nhalf+1):n)
    traceplot(as.mcmc(K[,i]), main=paste0("Observation ",i))

  title(paste0("Cluster allocation variables "), outer = TRUE)
  dev.off()
 
  # PSM
  K <- out$K
  psm <- PSM(K)
  png(file = paste0("PSM.png"), width = 8000, height = 5000, units = "px", res = 800)
  heatmap(psm, Rowv = NA, Colv = NA)
  dev.off()
  
  # estimate best partition
  part_BIN <- minbinder.ext(psm,cls.draw = K, method="all",include.greedy=TRUE)
  best.partition <- part_BIN$cl["best",]
  
  png(file = paste0("Partition.png"), width = 8000, height = 5000, units = "px", res = 800)
  matplot(t(X), type="l", col=best.partition)
  dev.off()
  
  # loglikelihood+counter
  logL      <- out$logL
  counter   <- out$counter
  png(file = paste0("LogL_counter.png"), width = 8000, height = 5000, units = "px", res = 800)
  par(mfrow=c(2,1))
  traceplot(as.mcmc(logL), main="Traceplot for the logLikelihood")
  traceplot(as.mcmc(counter), main="Traceplot of counter")
  text(350,15, labels=paste0('Overall proposed clusters: ',sum(counter)))
  dev.off()
  
  
  # print information on a txt file -----------------------------------------------------
  cat(
    "ALGORITHM\n",
    "Iter:   ", n_iter, "\n",
    "Burnin: ", burnin, "\n",
    "M:      ", M, "\n",
    "Mass:   ", mass,
    file="Output.txt", sep=""
  )
  
  setwd(dir.current)
}

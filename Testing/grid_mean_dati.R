setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
# setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)
library(roahd)
library(coda)

source("FBNP.R")
source("FBNP_hyper_alltime.R")
source("new_FBNP.R")
source("FBNP_hyper.R")
source("Prior Elicitation.R")
source('Smoothing.R')
source("new_Prior_elicitation.R")

#### DATA #### -------------------------------------------------------------------------------

# load data and rescale
load("X.RData")

matplot(t(X[12:16,]), type='l', lwd=1, lty=1, 
        main="", xlab="Time [ms]", #ylab="Evoked potential [micro volt]",
        ylab=TeX('Evoked Potential $\\[\\mu$V$\\]$'),
        #ylab=TeX('Evoked potential'),
        #ylab=paste0("Evoked potential [",expression(mu),"V]"),
        ylim=c(-700,650))

# eliminate bad data
#eliminate <- c(12,13,19,24)
eliminate <- c()
X <- X[-eliminate,] # matrix n x n_time, HO TOLTO ANCHE LA 24

# rescale data
#rescale <- 1 # 
rescale <- max(X)
X <- X/rescale 
matplot(t(X), type='l')

# cut x-axis
X_1 <- X[,seq(151,1050)]
matplot(t(X_1), type='l')
dim(X_1)
X <- X_1

smoothing_list <- smoothing(X = X, 
                            step = 9, 
                            nbasis = 25, 
                            spline_order = 4)

smoothing_list[['smoothing_parameters']][['rescale_parameter']] <- rescale
smoothing_list[['smoothing_parameters']][['observation_eliminated']] <- eliminate


###############################################################################################

setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
source("PSM.R")
library(coda)
library(devtools)
library(mcclust.ext)

#mean.phi.grid <- c(0.01,0.1,1,5,10,15,20)
mean.phi.grid <- c(1,2)
jj = 1

for(mean_phi in mean.phi.grid)
{
  print(jj)
  jj = jj+1
  
  hyper_list <- hyperparameters(var_phi = 0.01, 
                                X = smoothing_list$X,
                                beta = smoothing_list$beta,
                                scale = 1,
                                mean_phi = mean_phi)
  
  out <- FBNP_hyper(n_iter = 1000,
                    burnin = 0,
                    M = 1000,
                    mass = 10,
                    smoothing = smoothing_list,
                    hyperparam = hyper_list)
  
  run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                         'prior_parameters'     = hyper_list,
                         'smoothing_parameters' = smoothing_list$smoothing_parameters)
  
  
  # save old directory to get back there at the end
  dir.current <- getwd()
  
  # name of directory where I will put plots, I use current time in the name
  new.dir <- paste0(dir.current,"/Results/13Feb/FBNP/mass10_var_0.01_mean_phi_",mean_phi)
  
  # create such directory and go there
  dir.create(new.dir)
  setwd(new.dir) 
  
  # save RData ----------------------------------------------------------------------
  save(out, run_parameters, file = "Output.RData")
  
  # save traceplots -----------------------------------------------------------------
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
  
  # IG
  png(file = "InverseGamma.png", width = 8000, height = 5000, units = "px", res = 800)
  plot(seq(0,50,by=0.01), dinvgamma(seq(0,50,by=0.01),
                                    shape=hyper_list$c,
                                    rate=hyper_list$d)
  )
  dev.off()
  
  # faccio due plot del cluster assignment, così non mi da prolemi di margin too large
  nhalf <- round(n/2)
  
  # first half of observation
  png(file = "Traceplot_K_1.png", width = 8000, height = 5000, units = "px", res = 800)
  par(mfrow=n2mfrow(n/2))
  par(oma=c(0,0,2,0))
  for(i in 1:nhalf)
    traceplot(as.mcmc(K[-c(1:20),i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  #title("Cluster allocation variables", outer = TRUE)
  title(paste0("Cluster allocation variables "), outer = TRUE)
  dev.off()
  
  # second half of observations
  png(file = "Traceplot_K_2.png", width = 8000, height = 5000, units = "px", res = 800)
  par(mfrow=n2mfrow(n-nhalf))
  par(oma=c(0,0,2,0))
  for(i in (nhalf+1):n)
    traceplot(as.mcmc(K[-c(1:20),i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  #title("Cluster allocation variables", outer = TRUE)
  title(paste0("Cluster allocation variables "), outer = TRUE)
  dev.off()
  
  #for(i in 1:n)
  #{
  #  png(file = paste0("Obs_",i,".png"), width = 8000, height = 5000, units = "px", res = 800)
  #  traceplot(as.mcmc(K[,i]), main=paste0("Observation ",i))#, ylim=c(0,M))
  #  dev.off()
  #}
  
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
    "Mass:   ", mass, "\n\n",
    
    "SMOOTHING\n",
    "step:       ", run_parameters$smoothing_parameters$step, "\n",
    "L:          ", run_parameters$smoothing_parameters$number_basis, "\n",
    "rescale:    ", run_parameters$smoothing_parameters$rescale_parameter, "\n",
    "eliminated: ", run_parameters$smoothing_parameters$observation_eliminated,
    
    "var phi:    ", hyper_list$desired_prior_values,
    
    file="Output.txt", sep=""
  )
  
  setwd(dir.current)
  
}

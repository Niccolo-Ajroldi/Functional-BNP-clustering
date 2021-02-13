
# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
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

#### DATA ####-------------------------------------------------------------------------------

# simulate data from 2 Gaussian processes

n.1 <- 10
n.2 <- 10
n.3 <- 10
n <- n.1+n.2+n.3
n_time <- 100
time.grid <- seq(0, 10, length.out = n_time)

# Exponential covariance function over a time.grid

# tune correlation of simulated data:
# increase alpha to increase variability in each point
# increase beta to decrease covariance between times (high beta -> more rough function)
alpha <- 0.01
beta  <- 0.25
psi.1 <- exp_cov_function(time.grid, alpha, beta)

# mean function
mu.1 <- sin(0.2*pi*time.grid)
mu.2 <- sin(0.35*pi*(time.grid-4))
mu.3 <- sin(0.2*pi*(time.grid+2))

# simulate data
set.seed(1)
data.1 <- generate_gauss_fdata(n.1,mu.1,Cov=psi.1)
data.2 <- generate_gauss_fdata(n.2,mu.2,Cov=psi.1)
data.3 <- generate_gauss_fdata(n.3,mu.3,Cov=psi.1)

X <- rbind(data.1, data.2, data.3)
col <- c(rep(1,n.1), rep(2,n.2), rep(3,n.3))
matplot(time.grid, t(X), type='l', col=col, main="Simulated GP")

# rescale data
rescale <- 1 # rescale <- max(X)
X <- X/rescale 

# basis 
L <- 15
basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=L, norder=4)

# smooth data
X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)

# save coefficients
beta <- t(X_smoothed_f$fd$coefs)

# plot smoothed data
basis.t <- t(eval.basis(time.grid, basis))
X_smooth <- beta %*% basis.t
matplot(time.grid, t(X_smooth), type='l', col=col)

smoothing_parameters <- list('step' = 1,
                             'number_basis' = L,
                             'spline_order' = 4)
smoothing_list <- list('basis' = basis,
                       'beta'= beta,
                       'time.grid' = time.grid,
                       'X' = X,
                       'smoothing_parameters' = smoothing_parameters)

###############################################################################################

setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
source("PSM.R")
library(coda)
library(devtools)
library(mcclust.ext)

mass.grid <- c(0.5,1,5,10,50)
jj = 1

for(mass in mass.grid)
{
  print(jj)
  jj = jj+1
  
  hyper_list <- hyperparameters(var_phi = 1, 
                                X = smoothing_list$X,
                                beta = smoothing_list$beta,
                                scale = 1,
                                mean_phi = 10)
  
  out <- FBNP_hyper(n_iter = 200,
                    burnin = 0,
                    M = 500,
                    mass = mass,
                    smoothing = smoothing_list,
                    hyperparam = hyper_list)
  
  run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                         'prior_parameters'     = hyper_list,
                         'smoothing_parameters' = smoothing_list$smoothing_parameters)
  
  
  # save old directory to get back there at the end
  dir.current <- getwd()
  
  # name of directory where I will put plots, I use current time in the name
  new.dir <- paste0(dir.current,"/Results/TEST_13_2/FBNP/mean_10_var_1_mass_",mass)
  
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
                                    rate=hyper_list$d), type='l'
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

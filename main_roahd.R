setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
# setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)
library(latex2exp)
library(roahd)
source("FBNP.R")
source("FBNP_hyper.R")
source("FBNP_hyper_alltime.R")
source("Prior Elicitation.R")
source('Smoothing.R')
source("FBNP_orig_nosigma.R")

#### DATA #### -------------------------------------------------------------------------------

# load data and rescale
sani <- mfD_healthy
malati <- mfD_LBBB
x11()
plot(sani)
x11()
plot(malati)

sani.3 <- sani$fDList[[3]]$values
malati.3 <- malati$fDList[[3]]$values

X <- rbind(sani.3,malati.3)
idx <- c(rep(1,50),rep(2,50))

x11()
matplot(t(X), col=idx, type='l')

X <- X[,-seq(1,304)]

x11()
matplot(t(X), col=idx, type='l')

# eliminate bad data
#eliminate <- c(12,13,19,24)
eliminate <- c()
X <- X[-eliminate,] # matrix n x n_time, HO TOLTO ANCHE LA 24

# rescale data
#rescale <- 1 # 
rescale <- max(X)
X <- X/rescale
x11()
matplot(t(X), type='l', col=idx)

smoothing_list <- smoothing(X = X, 
                            step = 12, 
                            nbasis = 25, 
                            spline_order = 4)

smoothing_list[['smoothing_parameters']][['rescale_parameter']] <- rescale
smoothing_list[['smoothing_parameters']][['observation_eliminated']] <- eliminate



#### HYPERPARAM #### -------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(mean_phi = 10,
                              var_phi = 0.01, 
                              X = smoothing_list$X,
                              beta = smoothing_list$beta,
                              scale = 1)

#### CALL #### -------------------------------------------------------------------------------

out <- FBNP_hyper(n_iter = 500,
                  burnin = 0,
                  M = 500,
                  mass = 0.25,
                  smoothing = smoothing_list,
                  hyperparam = hyper_list)


### SAVE OUTPUT #### -------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters'     = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
)

out[['algorithm_parameters']] <- NULL # ok ma perché allora non salvarli direttamente da qui azichÃ¨ farli restiruire da FBNP e poi rimuoverli?

# save output
load('savez.R')
savez(out,'Big runs')


#### DIAGNOSTIC ####-------------------------------------------------------------------------

library(coda)
library(devtools)
library(mcclust.ext)

# traceplot of cluster allocation variables
source("traceplot_K.R")
traceplot_K(out, smoothing_list, run_parameters)  

x11()
traceplot(as.mcmc(out$counter))

source("PSM.R")
K <- out$K
psm <- PSM(K)
{x11(); heatmap(psm, Rowv = NA, Colv = NA)}
{x11(); plotpsm(psm)}



# estimate best partition
part_BIN <- minbinder.ext(psm,cls.draw = K, method="all",include.greedy=TRUE)
summary(part_BIN)
best.partition <- part_BIN$cl["best",]
x11()
matplot(t(X), type="l", col=best.partition)

save(out, X, file = "Results/Big runs/mass5_meanphi5_varphi0.01_M5000.RData")
load("Results/Big runs/mass5_meanphi5_varphi0.01_M5000.RData")

diagnostic(out, smoothing_list, run_parameters)
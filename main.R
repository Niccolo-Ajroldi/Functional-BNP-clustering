
# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
# setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)

# load FBNP function
source("FBNP.R")

# load function for prior elicitation
source("Prior Elicitation.R")



#### DATA #### -------------------------------------------------------------------------------

# load function smoothing
source('Smoothing.R')

# load data and rescale
load("X.RData")

#eliminates bad data
eliminate <- c(12,13,19,24)
X <- X[-eliminate,] # matrix n x n_time, HO TOLTO ANCHE LA 24


#rescale data
rescale <- 1
#rescale <- max(X)
X <- X/rescale 

matplot(t(X), type='l')

smoothing_list <- smoothing(X = X, 
                            step = 5, 
                            nbasis = 30, 
                            spline_order = 4)

smoothing_list[['smoothing_parameters']][['rescale_parameter']] <- rescale
smoothing_list[['smoothing_parameters']][['observation_eliminated']] <- eliminate

# save(X, basis, X_smoothed_f, beta, time.grid, file="Xdata.RData")

#### HYPERPARAM #### -------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(var_sigma = 100, var_phi = 100, 
                              X = smoothing_list$X,
                              beta = smoothing_list$beta)

# or set them a caso
# hyper_list <- list(a=2.1, b=1, c=2.1, d=1, m0=rep(0,L), Lambda0=diag(1,L))


#### CALL #### -------------------------------------------------------------------------------

out <- FBNP(n_iter = 1000,
            burnin = 500,
            thin = 1,
            M = 150,
            mass = 0.31,
            smoothing = smoothing_list,
            hyperparam = hyper_list)

### SAVE OUTPUT #### -------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters' = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
                       )

out[['algorithm_parameters']] <- NULL

# save output
save(out, run_parameters, file="Results/out_nico_24_12_eddajeee.RData")
save(out, run_parameters, file = "Results/out_10000iter_hyperacaso.RData")
save(out, run_parameters, file = "Ultimo.RData")





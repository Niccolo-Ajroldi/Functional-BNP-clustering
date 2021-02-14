# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)

source("FBNP.R")
source("Prior Elicitation.R")
source('Smoothing.R')

#### DATA #### -------------------------------------------------------------------------------

# load data and rescale
load("X.RData")

# eliminate bad data
eliminate <- c(12,13,19,24)
eliminate <- c()
X <- X[-eliminate,] # matrix n x n_time, HO TOLTO ANCHE LA 24

# rescale data
rescale <- 1 # rescale <- max(X)
X <- X/rescale 

matplot(t(X), type='l')

smoothing_list <- smoothing(X = X, 
                            step = 5, 
                            nbasis = 30, 
                            spline_order = 4)

smoothing_list[['smoothing_parameters']][['rescale_parameter']] <- rescale
smoothing_list[['smoothing_parameters']][['observation_eliminated']] <- eliminate

#### HYPERPARAM #### -------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(var_sigma = 10, var_phi = 10,
                              X = smoothing_list$X,
                              beta = smoothing_list$beta)

# or set them a caso
#L <- smoothing_list$smoothing_parameters$number_basis
#hyper_list <- list(a=2.1, b=1, c=2.1, d=1, m0=rep(0,L), Lambda0=diag(1,L))


#### CALL #### -------------------------------------------------------------------------------

out <- FBNP(n_iter = 300,
            burnin = 0,
            thin = 1,
            M = 150,
            mass = 0.31,
            smoothing = smoothing_list,
            hyperparam = hyper_list)

### SAVE OUTPUT #### -------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters'     = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
)

out[['algorithm_parameters']] <- NULL # ok ma perchè allora non salvarli direttamente da qui azichè farli restiruire da FBNP e poi rimuoverli?

# save output
#save(out, run_parameters, file = "Results/Aggiustato_2_01.RData")

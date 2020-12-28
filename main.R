
#setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
#setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)

# load FBNP function
source("FBNP.R")

# load function for prior elicitation
source("Prior Elicitation.R")

#### DATA #### -------------------------------------------------------------------------------

step      <- 5
n_time    <- 1600/step
time.grid <- 1:n_time # TODO: riscalare la time.grid
time.grid <- time.grid#/max(time.grid)

L <- 30
m <- 4

# load data and rescale
load("Smoothing/smooth_60b_nopenalization.RData")
X <- X[-c(12,13,19,24),] # matrix n x n_time, HO TOLTO ANCHE LA 24
X <- X[, seq(1,1600,by=step)]
X <- 10*X/max(X) # TODO: capire se si può migliorare

matplot(t(X), type='l')

# basis 
basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=L, norder=m)

# smoothing
X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)

# save coefficients
beta <- t(X_smoothed_f$fd$coefs)

#save(X, basis, X_smoothed_f, beta, time.grid, file="Xdata.RData")

#### HYPERPARAM #### -------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(var_sigma=10, var_phi=10, X, beta)

# or set them a caso
#hyper_list <- list(a=2.1, b=1, c=2.1, d=1, m0=rep(0,L), Lambda0=diag(1,L))


#### CALL #### -------------------------------------------------------------------------------

out <- FBNP(n_iter = 5000,
            burnin = 3000,
            thin = 1,
            M = 150,
            mass = 0.31,
            X = X,
            basis = basis,
            beta = beta,
            time.grid = time.grid,
            hyperparam = hyper_list)


# save output
save(out, file="Results/out_nico_24_12_eddajeee.RData")
save(out, file = "Results/out_10000iter_hyperacaso.RData")





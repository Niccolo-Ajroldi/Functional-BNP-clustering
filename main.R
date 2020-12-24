
#setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
#setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)

# load FBNP function
source("FBNP.R")

#### DATA #### -------------------------------------------------------------------------------

step      <- 10
n_time    <- 1600/step
time.grid <- 1:(1600/step) # TODO: riscalare la time.grid

L <- 30
m <- 4

# load data and rescale
load("Smoothing/smooth_60b_nopenalization.RData")
X <- X[-c(12,13,19),] # matrix n x n_time
X <- X[, seq(1,1600,by=step)]
X <- X/max(X) # TODO: capire se si può migliorare

# basis 
basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=L, norder=m)

# smoothing
X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)

# save coefficients
beta <- t(X_smoothed_f$fd$coefs)

#### CALL #### -------------------------------------------------------------------------------

out <- FBNP(n_iter = 1000,
            burnin = 500,
            thin = 1,
            M = 150,
            mass = 100,
            X = X,
            basis = basis,
            beta = beta,
            time.grid = time.grid)


#### DIAGNOSTIC #### -------------------------------------------------------------------------

# save output
# save(out, file="Results/out_nico_23_12_500_iter_.RData") 

names(out)

K            <- out$K
mu_coef_out  <- out$mu_coef_out
sigma2_out   <- out$sigma2_out
probs_j_out  <- out$probs_j_out
probs_ij_out <- out$probs_ij_out


# perform diagnostic
source("diagnostic.R")



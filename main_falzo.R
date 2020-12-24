
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

L <- 30
m <- 4

# load data
data(kma.data)
X <- kma.data$y0 # matrix n x n_time
time.grid <- (kma.data$x)[1,] # time grid
n_time <- length(time.grid)

# basis
basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=L, norder=m)

# smoothing
X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)

# save coefficients
beta <- t(X_smoothed_f$fd$coefs)

#### CALL #### --------------------------------------------------------------------------

K <- FBNP(n_iter = 2000,
          burnin = 1000,
          thin = 1,
          M = 150,
          mass = 10,
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









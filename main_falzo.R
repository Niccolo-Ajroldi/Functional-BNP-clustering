
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

#### HYPERPARAM #### -------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(var_sigma=10, var_phi=10, X, beta)

# or set them a caso
#hyper_list <- list(a=2.1, b=1, c=2.1, d=1, m0=rep(0,L), Lambda0=diag(1,L))


#### CALL #### --------------------------------------------------------------------------

out <- FBNP(n_iter = 500,
            burnin = 100,
            thin = 1,
            M = 150,
            mass = 1,
            X = X,
            basis = basis,
            beta = beta,
            time.grid = time.grid,
            hyperparam = hyper_list)


#### DIAGNOSTIC #### -------------------------------------------------------------------------

# save output
save(out, file="Results/out_nico_falZo_27_12_2000_iter.RData") 

names(out)

K            <- out$K
mu_coef_out  <- out$mu_coef_out
sigma2_out   <- out$sigma2_out
probs_j_out  <- out$probs_j_out
probs_ij_out <- out$probs_ij_out


# perform diagnostic
source("diagnostic.R")









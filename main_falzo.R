
#setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
#setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)

source("FBNP.R")
source("Prior Elicitation.R")
source('Smoothing.R')

#### DATA #### -------------------------------------------------------------------------------

# load data
data(kma.data)
X <- kma.data$y0 # matrix n x n_time
time.grid <- (kma.data$x)[1,] # time grid
#time.grid <- 1:200

n_time <- length(time.grid)

matplot(t(X), type='l')

# basis 
L <- 30
basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=L, norder=4)

# smooth data
X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)

# save coefficients
beta <- t(X_smoothed_f$fd$coefs)

smoothing_parameters <- list('step' = 1,
                             'number_basis' = L,
                             'spline_order' = 4)

smoothing_list <- list('basis' = basis,
                       'beta'= beta,
                       'time.grid' = time.grid,
                       'X' = X,
                       'smoothing_parameters' = smoothing_parameters)

#### HYPERPARAM #### -------------------------------------------------------------------------------

# elicit hyperparameters
#hyper_list <- hyperparameters(var_sigma = 100, var_phi = 100, 
#                              X = smoothing_list$X,
#                              beta = smoothing_list$beta)
#
# or set them a caso
hyper_list <- list(a=2.01, b=1, c=2.01, d=1, m0=rep(0,L), Lambda0=diag(1,L))


#### CALL #### --------------------------------------------------------------------------

out <- FBNP(n_iter = 3000,
            burnin = 2000,
            thin = 1,
            M = 150,
            mass = 1000,
            smoothing = smoothing_list,
            hyperparam = hyper_list)



### SAVE OUTPUT #### -------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters' = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
                       )

out[['algorithm_parameters']] <- NULL


#### DIAGNOSTIC #### -------------------------------------------------------------------------

# save output
#save(out, file="Results/out_nico_falZo_27_12_2000_iter.RData") 

names(out)

K            <- out$K
mu_coef_out  <- out$mu_coef_out
sigma2_out   <- out$sigma2_out
probs_j_out  <- out$probs_j_out
probs_ij_out <- out$probs_ij_out


K            <- out$K
View(K)








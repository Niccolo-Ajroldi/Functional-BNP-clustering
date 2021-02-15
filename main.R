
# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
#setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
 setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

rm(list=ls())
cat("\014")

library(fda)
library(latex2exp)
source("FBNP.R")
source("FBNP_hyper.R")
source("hyperparameters.R")
source('Smoothing.R')

#### DATA #### -------------------------------------------------------------------------------

# load data and rescale
load("X.RData")

matplot(t(X[12:16,]), type='l', lwd=1, lty=1, 
        main="", xlab="Time [ms]",
        ylab=TeX('Evoked Potential $\\[\\mu$V$\\]$'),
        ylim=c(-700,650))

# rescale data
rescale <- max(X)
X <- X/rescale 

# cut x-axis
X_1 <- X[,seq(151,1050)]
matplot(t(X_1), type='l')
dim(X_1)
X <- X_1

smoothing_list <- smoothing(X = X, 
                            step = 12, 
                            nbasis = 20, 
                            spline_order = 4)


#### HYPERPARAM #### -------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(X = smoothing_list$X,
                              beta = smoothing_list$beta,
                              mean_phi = 500,
                              var_phi = 10)

#### CALL #### -------------------------------------------------------------------------------

out <- FBNP_hyper(n_iter = 35000,
                  burnin = 10000,
                  M = 500,
                  mass = 0.5,
                  smoothing = smoothing_list,
                  hyperparam = hyper_list)


### SAVE OUTPUT #### -------------------------------------------------------------------------

save(out, smoothing_list, hyper_list, file = "Results/Big runs/mass100_meanphi0.9_varphi0.001_M1000.RData")

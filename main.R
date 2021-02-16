rm(list=ls())
cat("\014")

library(fda)
library(latex2exp)
source("Tools/FBNP/FBNP.R")
source("Tools/FBNP/FBNP_hyper.R")
source("Tools/hyperparameters.R")
source('Tools/Smoothing.R')

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

out <- FBNP_hyper(n_iter = 5,
                  burnin = 0,
                  M = 500,
                  mass = 0.5,
                  smoothing = smoothing_list,
                  hyperparam = hyper_list)


### SAVE OUTPUT #### -------------------------------------------------------------------------

save(out, smoothing_list, hyper_list, file = "Last_run.RData")

source("Tools/save_fun.R")
save_fun(out,'Last_run')

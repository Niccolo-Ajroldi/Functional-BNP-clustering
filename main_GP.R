
# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)
library(roahd)

source("FBNP.R")
source("FBNP_hyper_alltime.R")
source("new_FBNP.R")
source("FBNP_hyper.R")
source("Prior Elicitation.R")
source('Smoothing.R')
source("new_Prior_elicitation.R")

#### DATA ####-------------------------------------------------------------------------------

# simulate data from 2 Gaussian processes

n.1 <- 15
n.2 <- 15
n <- n.1+n.2
#n_time <- 300
time.grid <- seq(0, 10, length.out = 50)

# Exponential covariance function over a time.grid

# tune correlation of simulated data:
# increase alpha to increase variability in each point
# increase beta to decrease covariance between times (high beta -> more rough function)
alpha <- 0.1
beta  <- 0.5
psi.1 <- exp_cov_function(time.grid, alpha, beta)
psi.2 <- psi.1
image(psi.1,
      main = 'Exponential covariance function',
      xlab = 'time.grid', ylab = 'time.grid')

# mean function
mu.1 <- sin(0.2*pi*time.grid)
mu.2 <- sin(0.35*pi*(time.grid-4))
plot(time.grid, mu.1, ylab = "y(t)", type='l', col=1)
lines(time.grid, mu.2, ylab = "y(t)", type='l', col=2)

#to simulate data you use...
set.seed(1)
data.1 <- generate_gauss_fdata(n.1,mu.1,Cov=psi.1)
data.2 <- generate_gauss_fdata(n.2,mu.2,Cov=psi.1)

X <- rbind(data.1, data.2)
col <- c(rep(1,n.1), rep(2,n.2))
matplot(time.grid, t(X), type='l', col=col, main="Simulated GP")

# rescale data
rescale <- 1 # rescale <- max(X)
X <- X/rescale 

# basis 
L <- 15
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

#### HYPERPARAM ####-------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(var_phi = 1e5, 
                              X = smoothing_list$X,
                              beta = smoothing_list$beta,
                              scale = 1)


#### CALL ####--------------------------------------------------------------------------

out <- FBNP_hyper(n_iter = 3000,
                          burnin = 0,
                          M = 2000,
                          mass = 0.5,
                          smoothing = smoothing_list,
                          hyperparam = hyper_list)


### SAVE OUTPUT ####-------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters' = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
)

out[['algorithm_parameters']] <- NULL

save(out, file="Results/out_post_corrado.RData") 

#### DIAGNOSTIC ####-------------------------------------------------------------------------

library(coda)
library(devtools)
library(mcclust.ext)

# traceplot of cluster allocation variables
source("traceplot_K.R")
traceplot_K(out, smoothing_list, run_parameters)  

n.1+n.2
apply(out$K[-1,], 2, max)


# posterior similarity matrix
source("PSM.R")
K <- out$K
psm <- PSM(K)
{x11(); heatmap(psm, Rowv = NA, Colv = NA)}
{x11(); plotpsm(psm)}

# estimate best partition
part_BIN <- minbinder.ext(psm,cls.draw = K, method="all",include.greedy=TRUE)
summary(part_BIN)

# Plot of the partitions with the different methods
x11()
par(mfrow=c(2,3))
for(ii in 1:5)
  matplot(t(X), col=(part_BIN$cl)[ii,], type='l', main=(row.names(part_BIN$cl))[ii])


# choose a single partition, in this case "avg"
partition.BIN <- minbinder.ext(psm,cls.draw = K, method="draws")[[1]]

x11()
matplot(t(X), type="l", col=partition.BIN)



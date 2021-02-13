
setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
# setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
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
source("FBNP_orig_nosigma.R")

#### DATA ####-------------------------------------------------------------------------------

# simulate data from 2 Gaussian processes

n.1 <- 10
n.2 <- 10
n.3 <- 10
n <- n.1+n.2+n.3
n_time <- 100
time.grid <- seq(0, 10, length.out = n_time)

# Exponential covariance function over a time.grid

# tune correlation of simulated data:
# increase alpha to increase variability in each point
# increase beta to decrease covariance between times (high beta -> more rough function)
alpha <- 0.05
beta  <- 0.5
psi.1 <- exp_cov_function(time.grid, alpha, beta)
psi.3 <- psi.2 <- psi.1
image(psi.1,
      main = 'Exponential covariance function',
      xlab = 'time.grid', ylab = 'time.grid')

# mean function
mu.1 <- sin(0.2*pi*time.grid)
mu.2 <- sin(0.35*pi*(time.grid-4))
mu.3 <- sin(0.25*pi*(time.grid-2))
plot(time.grid, mu.1, ylab = "y(t)", type='l', col=1)
lines(time.grid, mu.2, ylab = "y(t)", type='l', col=2)
lines(time.grid, mu.3, ylab = "y(t)", type='l', col=3)

# simulate data
set.seed(1)
data.1 <- generate_gauss_fdata(n.1,mu.1,Cov=psi.1)
data.2 <- generate_gauss_fdata(n.2,mu.2,Cov=psi.1)
data.3 <- generate_gauss_fdata(n.3,mu.3,Cov=psi.1)

X <- rbind(data.1,data.2)
col <- c(rep(1,n.1), rep(2,n.2))
#X <- rbind(data.1, data.2, data.3)
#col <- c(rep(1,n.1), rep(2,n.2), rep(3,n.3))
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

# plot smoothed data
basis.t <- t(eval.basis(time.grid, basis))
X_smooth <- beta %*% basis.t
matplot(time.grid, t(X_smooth), type='l', col=col)

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
hyper_list <- hyperparameters(mean_phi=10,
                              var_phi = 0.01, 
                              X = smoothing_list$X,
                              beta = smoothing_list$beta,
                              scale = 100)


#### CALL ####--------------------------------------------------------------------------

out <- FBNP_hyper(n_iter = 100,
                          burnin = 0,
                          M = 1000,
                          mass = 0.5,
                          smoothing = smoothing_list,
                          hyperparam = hyper_list)


### SAVE OUTPUT ####-------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters' = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
)

out[['algorithm_parameters']] <- NULL

#save(out, file="Results/nico_11_2")
save(out, file="Results/tere_orig_nosigma_m100v1e3")

#### DIAGNOSTIC ####-------------------------------------------------------------------------

library(coda)
library(devtools)
library(mcclust.ext)
load("Results/tere_orig_nosigma_3gr_v1e-1")
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
best.partition <- part_BIN$cl["best",]
x11()
matplot(t(X), type="l", col=best.partition)

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



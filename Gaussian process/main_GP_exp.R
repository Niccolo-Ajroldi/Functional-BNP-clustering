rm(list=ls())
cat("\014")

library(fda)
library(fdakma)
library(roahd)
library(coda)

source("Tools/FBNP/FBNP.R")
source("Tools/FBNP/FBNP_hyper.R")
source("Tools/hyperparameters.R")
source('Tools/Smoothing.R')

#### DATA ####-------------------------------------------------------------------------------

# simulate data from 3 Gaussian processes
n.1 <- 10
n.2 <- 10
n.3 <- 10
n <- n.1+n.2+n.3

# time grid
n_time <- 100
time.grid <- seq(0, 10, length.out = n_time)

library(invgamma)

# diagonal covariance matrix
mean_phi=7
var_phi=4
c <- mean_phi^2/var_phi + 2
d <- mean_phi^3/var_phi + mean_phi
psi.1 <- diag(nrow=n_time,ncol=n_time) * rinvgamma(n_time, shape=c, rate=d)
image(psi.1)

# mean function
mu.1 <- 10*sin(0.2*pi*time.grid)
mu.2 <- 10*sin(0.35*pi*(time.grid-4))
mu.3 <- 10*sin(0.2*pi*(time.grid+2))

# bind simulated data
set.seed(1)
data.1 <- generate_gauss_fdata(n.1,mu.1,Cov=psi.1)
data.2 <- generate_gauss_fdata(n.2,mu.2,Cov=psi.1)
data.3 <- generate_gauss_fdata(n.3,mu.3,Cov=psi.1)
X <- rbind(data.1, data.2, data.3)

# plot simulated data
col <- c(rep(1,n.1), rep(2,n.2), rep(3,n.3))
matplot(time.grid, t(X), type='l', col=col, main="Simulated GP")

#png(file = paste0("SimulatedGP.png"), width = 8000, height = 5000, units = "px", res = 800)
#matplot(time.grid, t(X), type='l', col=col, main="Simulated GP", ylab="y")
#dev.off()

# rescale data
rescale <- max(X)
X <- X/rescale 

# number of data
n <- dim(X)[1]

# Fourier basis
basis <- create.fourier.basis(rangeval=range(time.grid), nbasis=7)

# coefficients
beta <- t(Data2fd(y=t(X), argvals=time.grid, basisobj=basis)$coefs)

# smooth data
basis.t <- t(eval.basis(time.grid, basis))
X_smooth <- beta %*% basis.t

# plot smoothed data
matplot(time.grid, t(X_smooth), type='l', col=col)

smoothing_parameters <- list('step' = 1,
                             'number_basis' = 7,
                             'spline_order' = 1)
smoothing_list <- list('basis' = basis,
                       'beta'= beta,
                       'time.grid' = time.grid,
                       'X' = X,
                       'smoothing_parameters' = smoothing_parameters)

#### HYPERPARAM ####-------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(var_phi = 0.001, 
                              X = smoothing_list$X,
                              beta = smoothing_list$beta,
                              mean_phi = 1.1)


#### CALL ####-------------------------------------------------------------------------------

out <- FBNP_hyper(n_iter = 8000,
                          burnin = 3000,
                          M = 500,
                          mass = 65,
                          smoothing = smoothing_list,
                          hyperparam = hyper_list)


### SAVE OUTPUT ####-------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters' = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
)

out[['algorithm_parameters']] <- NULL

source("savez.R")
savez(out, "/Results/GP_exp_def")

#save(out, file="Results/nico_11_2")
#save(out, file="Results/tere_orig_nosigma_m100v1e3")

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



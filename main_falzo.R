
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

png(file = "Falzo.png", width = 4500, height = 4000, units = "px", res = 1000)
matplot(time.grid, t(X), type='l', lwd=1, lty=2, 
        main="", xlab="x", ylab="y")
dev.off()

n_time <- length(time.grid)

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
hyper_list <- hyperparameters(var_sigma = 10, var_phi = 10, 
                              X = smoothing_list$X,
                              beta = smoothing_list$beta)

# or set them a caso
#hyper_list <- list(a=2.1, b=1, c=2.1, d=1, m0=rep(0,L), Lambda0=diag(1,L))


#### CALL #### --------------------------------------------------------------------------

out <- FBNP(n_iter = 5000,
            burnin = 3000,
            thin = 1,
            M = 150,
            mass = 0.7, # voglio provare con la massa alta
            smoothing = smoothing_list,
            hyperparam = hyper_list)



### SAVE OUTPUT #### -------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters' = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
                       )

out[['algorithm_parameters']] <- NULL

save(out, file="Results/out_nico_falZo_2_1.RData") 


#### DIAGNOSTIC #### -------------------------------------------------------------------------


library(coda)
library(devtools)
library(mcclust.ext)

K <- out$K

source("PSM.R")
psm <- PSM(K)

x11()
heatmap(psm, Rowv = NA, Colv = NA)

est_part_BINDER <- function(clust, PSM){
  res_loss <- c()
  for(i in 1:nrow(clust)){
    # rappresentazione binaria della partizione
    binary_part <- sapply(clust[i,], function(x) as.numeric(x == clust[i,]))
    # loss per la partizione i
    res_loss[i] <- sum(abs(binary_part - PSM))
  }
  # tra le iterazioni trova quella che minimizza
  return(clust[which.min(res_loss),])
}

part_BIN <- as.numeric(as.factor( est_part_BINDER(K,PSM(K)) ))
table(part_BIN)

matplot(t(X), type="l", col=part_BIN)


library(coda)
library(devtools)
library(mcclust.ext)

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





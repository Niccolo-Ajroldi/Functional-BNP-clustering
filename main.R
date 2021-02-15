
# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
#setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
 setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

rm(list=ls())
cat("\014")

library(fda)
library(fdakma)
library(latex2exp)
source("FBNP.R")
source("FBNP_hyper.R")
source("FBNP_hyper_alltime.R")
source("hyperprior.R")
source('Smoothing.R')
source("FBNP_orig_nosigma.R")

#### DATA #### -------------------------------------------------------------------------------

# load data and rescale
load("X.RData")

matplot(t(X[12:16,]), type='l', lwd=1, lty=1, 
        main="", xlab="Time [ms]", #ylab="Evoked potential [micro volt]",
        ylab=TeX('Evoked Potential $\\[\\mu$V$\\]$'),
        #ylab=TeX('Evoked potential'),
        #ylab=paste0("Evoked potential [",expression(mu),"V]"),
        ylim=c(-700,650))

# eliminate bad data
#eliminate <- c(12,13,19,24)
eliminate <- c()
X <- X[-eliminate,] # matrix n x n_time, HO TOLTO ANCHE LA 24

# rescale data
#rescale <- 1 # 
#rescale <- max(X)
#X <- X/rescale 
#matplot(t(X), type='l')

png(file = paste0("Data_raw.png"), width = 8000, height = 5000, units = "px", res = 800)
matplot(t(X), type='l', lwd=1, lty=1, 
        xlab="Time [ms]",
        ylab=TeX('Evoked Potential $\\[\\mu$V$\\]$'),
        main="Evoked potential")
dev.off()

# cut x-axis
X_1 <- X[,seq(151,1050)]
matplot(t(X_1), type='l')
dim(X_1)
X <- X_1

png(file = paste0("Data_cutted.png"), width = 8000, height = 5000, units = "px", res = 800)
matplot(seq(151,1050),
        t(X), type='l', lwd=1, lty=1, 
        xlab="Time [ms]",
        ylab=TeX('Evoked Potential $\\[\\mu$V$\\]$'),
        main="Evoked potential")
dev.off()

smoothing_list <- smoothing(X = X, 
                            step = 12, 
                            nbasis = 20, 
                            spline_order = 4)

smoothing_list[['smoothing_parameters']][['rescale_parameter']] <- rescale
smoothing_list[['smoothing_parameters']][['observation_eliminated']] <- eliminate

#### HYPERPARAM #### -------------------------------------------------------------------------------

# elicit hyperparameters
hyper_list <- hyperparameters(mean_phi = 0.85,
                              var_phi = 0.001, 
                              X_obs = smoothing_list$X,
                              beta = smoothing_list$beta)

#### CALL #### -------------------------------------------------------------------------------

out <- FBNP_hyper(n_iter = 8000,
                  burnin = 3000,
                  M = 500,
                  mass = 60,
                  smoothing = smoothing_list,
                  hyperparam = hyper_list)


### SAVE OUTPUT #### -------------------------------------------------------------------------

run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                       'prior_parameters'     = hyper_list,
                       'smoothing_parameters' = smoothing_list$smoothing_parameters
                       )

out[['algorithm_parameters']] <- NULL # ok ma perché allora non salvarli direttamente da qui azichÃ¨ farli restiruire da FBNP e poi rimuoverli?

source('savez.R')
savez(out,'TEST_15_2/mass60_mean085_rifatto')

#### DIAGNOSTIC ####-------------------------------------------------------------------------

#rm(list=ls())
#setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
#load("Results/NicoNight/mass60_mean085/Output.RData")

library(coda)
library(devtools)
library(mcclust.ext)

# traceplot of cluster allocation variables
source("traceplot_K.R")
traceplot_K(out, smoothing_list, run_parameters)  

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

#save(out, X, file = "Results/Big runs/mass5_meanphi5_varphi0.01_M5000.RData")
#load("Results/Big runs/mass5_meanphi5_varphi0.01_M5000.RData")

diagnostic(out, smoothing_list, run_parameters)

# Plot of the partitions with the different methods
x11()
par(mfrow=c(2,3))
for(ii in 1:5)
        matplot(t(X), col=(part_BIN$cl)[ii,], type='l', main=(row.names(part_BIN$cl))[ii])


# choose a single partition, in this case "avg"
partition.BIN <- minbinder.ext(psm,cls.draw = K, method="comp")[[1]]
x11()
matplot(t(X), type="l", col=partition.BIN)


load("GOSE.RData")
GOSE.rm <- GOSE

partition <- partition.BIN

# partition by clusters
x11()
par(mfrow = c(3,3))
count <- 1
for(j in levels(as.factor(partition)))
{
        indexes.j <- which(partition==j)
        if(length(indexes.j)==1)
                matplot(X[indexes.j,], type='l', col=GOSE.rm[indexes.j], xlab="", main=paste0("Cluster ",count))
        else
                matplot(t(X[indexes.j,]), type='l', col=GOSE.rm[indexes.j], xlab="", main=paste0("Cluster ",count))
        count <- count+1
}

# partition by GOSE
x11()
par(mfrow = c(2,1))
for(i in 1:2)
{
        matplot(t(X[GOSE.rm==i,]), col=partition[GOSE.rm==i], type = 'l', main = paste("GOSE",i), xlab="")
}



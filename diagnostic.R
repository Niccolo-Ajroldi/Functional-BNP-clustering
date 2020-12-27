setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")

rm(list=ls())
cat("\014")

#### DIAGNOSTIC #####------------------------------------------

library(coda)
library(devtools)
# devtools::install_github("sarawade/mcclust.ext")
library(mcclust.ext)
library(fda)
library(fdakma)


load("Results/out_10000iter_hyperacaso.RData")
load("Xdata.RData")
load("GOSEcomparison.RData")

#### save in list observations with GOSE=1 and GOSE=2 ####
observations.GOSE <- list()
for(i in 1:2)
{
  group <- X[which(GOSE.rm==i),]
  observations.GOSE[[i]] <- group
}

#### output of the algorithm ####
names(out)

K            <- out$K
mu_coef_out  <- out$mu_coef_out
sigma2_out   <- out$sigma2_out
probs_j_out  <- out$probs_j_out
probs_ij_out <- out$probs_ij_out

# Evaluation of the Posterior Similarity Matrix
source("PSM.R")
psm <- PSM(K)

#### Number of clusters ####

unique_clusters <- apply(K, 1, function(x) length(unique(x)))
coda_import <- as.mcmc(unique_clusters)
summary(coda_import)
geweke.diag(coda_import)
acfplot(coda_import)

x11()
heatmap(psm, Rowv = NA, Colv = NA)

x11()
plotpsm(psm)
# non capisco perché vengano diverse, dato che entrambe dovrebbero funzionare
# usando default values di dist (distance = "eucledian") e hclust (method = "complete")

#### Partition BINDER - using mcclust.ext ####

part_BIN <- minbinder.ext(psm,cls.draw = K, method="all",include.greedy=TRUE)
summary(part_BIN)

partition.BIN <- minbinder.ext(psm,cls.draw = K, method="avg")[[1]]

# CREDIBLE BALL
# The credible ball summarizes the uncertainty in the posterior around a
# clustering estimate part_VI and is defined as the smallest
# ball around part_VI with posterior probability at least 1-alpha

credibleball.BIN <- credibleball(partition.BIN, K, c.dist = "BIN", alpha = 0.05)
summary(credibleball.BIN)


#### Partition VI - using mcclust.ext####

part_VI <- minVI(psm, cls.draw=K, method=("all"), include.greedy=TRUE,
                 start.cl=NULL, suppress.comment=TRUE)
summary(part_VI)

partition.VI <- minVI(psm, cls.draw=K, method=("avg"),suppress.comment=TRUE)[[1]]

credibleball.VI <- credibleball(partition.VI, K, c.dist = "VI", alpha = 0.05)
summary(credibleball.VI)


#### PLOT ####------------------------------------------------------------

# Since the best partition with VI pushes all observations in the same cluster,
# here I do a little analysis with the results of the Binder method
# Plots of the best partitions of the partitions in the credible ball

## BEST PARTITION
partition <- as.numeric(credibleball.BIN$c.star)

x11()
matplot(t(X), col=partition, type='l')

# save observations per clusters
indexes <- unique(partition)
observations <- list()
idx=1
for(i in indexes)
{
  clust <- X[which(partition==i),]
  observations[[idx]] <- clust
  idx=idx+1
}

table(partition)

x11()
par(mfrow = c(3,2))
for(i in 1:idx)
{
  if(length(observations[[i]]) > length(time.grid))
  {
    matplot(t(observations[[i]]), type = 'l', col=GOSE.rm, main = paste("Cluster",i), xlab="")
  }else{
    matplot(observations[[i]], type = 'l', col=GOSE.rm, main = paste("Cluster",i), xlab="")
    }
}

# partition by GOSE
x11()
par(mfrow = c(2,1))
for(i in 1:2)
{
  matplot(t(observations.GOSE[[i]]), type = 'l', main = paste("GOSE",i), xlab="")
}

# We may repeat the plot for other partitions in the credible ball
partition <- as.numeric(credibleball.BIN$c.lowervert)
partition <- as.numeric(credibleball.BIN$c.horiz)


#### plot mu ####--------------------------------------------------



  

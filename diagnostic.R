
# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

cat("\014")

library(coda)
library(devtools)
library(mcclust.ext)

#### DATA #####----------------------------------------------------------------------------

# load results
#load("Results/out_10000iter_hyperacaso.RData")

# load X
load("X.RData")

# load complete GOSE
load("GOSE.RData")

# deleted observations
eliminated <- run_parameters$smoothing_parameters$observation_eliminated

# remove deleted observations from gose
GOSE.rm <- GOSE[-eliminate]

# output of the algorithm
K            <- out$K
mu_coef_out  <- out$mu_coef_out
sigma2_out   <- out$sigma2_out
probs_j_out  <- out$probs_j_out
probs_ij_out <- out$probs_ij_out

#### PSM #####----------------------------------------------------------------------------

# Evaluation of the Posterior Similarity Matrix
source("PSM.R")
psm <- PSM(K)

x11()
heatmap(psm, Rowv = NA, Colv = NA)

x11()
plotpsm(psm)

# Number of clusters
unique_clusters <- apply(K, 1, function(x) length(unique(x)))
coda_import <- as.mcmc(unique_clusters)
summary(coda_import)
geweke.diag(coda_import)
acfplot(coda_import)


#### Partition BINDER - using mcclust.ext ####--------------------------------------------

part_BIN <- minbinder.ext(psm,cls.draw = K, method="all",include.greedy=TRUE)
summary(part_BIN)

# Plot of the partitions with the different methods
x11()
par(mfrow=c(2,3))
for(ii in 1:5)
  matplot(t(X), col=(part_BIN$cl)[ii,], type='l', main=(row.names(part_BIN$cl))[ii])


# choose a single partition, in this case "avg"
partition.BIN <- minbinder.ext(psm,cls.draw = K, method="avg")[[1]]

# CREDIBLE BALL
# The credible ball summarizes the uncertainty in the posterior around a
# clustering estimate part_VI and is defined as the smallest
# ball around part_VI with posterior probability at least 1-alpha

credibleball.BIN <- credibleball(partition.BIN, K, c.dist = "BIN", alpha = 0.05)
summary(credibleball.BIN)


#### Partition VI - using mcclust.ext####-----------------------------------------------------

part_VI <- minVI(psm, cls.draw=K, method=("all"), include.greedy=TRUE,
                 start.cl=NULL, suppress.comment=TRUE)
summary(part_VI)

partition.VI <- minVI(psm, cls.draw=K, method=("avg"),suppress.comment=TRUE)[[1]]

credibleball.VI <- credibleball(partition.VI, K, c.dist = "VI", alpha = 0.05)
summary(credibleball.VI)


#### PLOT of Best partition ####------------------------------------------------------------

# partition to plot
partition <- partition.BIN
length(unique(partition))

x11()
matplot(t(X), col=partition, type='l')

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


#### PLOT of partitions in credible balls ####--------------------------------------------------

# Since the best partition with VI pushes all observations in the same cluster,
# here I do a little analysis with the results of the Binder method
# Plots of the best partitions of the partitions in the credible ball

partition <- as.numeric(credibleball.BIN$c.star)

# partition by clusters
x11()
par(mfrow = c(3,2))
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

# We may repeat the plot for other partitions in the credible ball
#partition <- as.numeric(credibleball.BIN$c.lowervert)
#partition <- as.numeric(credibleball.BIN$c.horiz)



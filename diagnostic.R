setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
# setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')


library(coda)
library(devtools)
library(mcclust.ext)


#### DIAGNOSTIC ####-------------------------------------------------------------------------

load("Results/Big runs/mass100_meanphi0.9_varphi0.001_M1000.RData")

# traceplot of cluster allocation variables
source("traceplot_K.R")
traceplot_K(out, smoothing_list)  

source("PSM.R")
K <- out$K
psm <- PSM(K)
{x11(); heatmap(psm, Rowv = NA, Colv = NA)}
{x11(); plotpsm(psm)}


#### PARTITION ANALYSIS ####----------------------------------------------------------------

## Binder Loss Function
part_BIN <- minbinder.ext(psm,cls.draw = K, method="all",include.greedy=TRUE)
summary(part_BIN)

# Plot of the partitions with the different methods
x11()
par(mfrow=c(2,3))
for(ii in 1:5)
  matplot(t(X), col=(part_BIN$cl)[ii,], type='l', main=(row.names(part_BIN$cl))[ii])

# Plot the best partition
best.partition <- part_BIN$cl["best",]
x11()
matplot(t(X), type="l", col=best.partition)


## Variation of Information
part_VI <- minVI(psm, cls.draw=K, method=("all"), include.greedy=TRUE,
                 start.cl=NULL, suppress.comment=TRUE)
summary(part_VI)

# Plot of the partitions with the different methods
x11()
par(mfrow=c(2,3))
for(ii in 1:5)
  matplot(t(X), col=(part_VI$cl)[ii,], type='l', main=(row.names(part_VI$cl))[ii])

# Plot the best partition
best.partition <- part_VI$cl["best",]
x11()
matplot(t(X), type="l", col=best.partition)


#### FACTOR comparison #####----------------------------------------------------------------------------

# load complete GOSE
load("GOSE.RData")
load("LCF.RData")
load("DRS.RData")

partition <- best.partition

# factor to investigate
fac <- GOSE

# partition by clusters
x11()
par(mfrow = c(3,3))
count <- 1
for(j in levels(as.factor(partition)))
{
  indexes.j <- which(partition==j)
  if(length(indexes.j)==1)
    matplot(X[indexes.j,], type='l', col=fac[indexes.j], xlab="", main=paste0("Cluster ",count))
  else
    matplot(t(X[indexes.j,]), type='l', col=fac[indexes.j], xlab="", main=paste0("Cluster ",count))
  count <- count+1
}

# partition by factor
x11()
par(mfrow = c(2,1))
for(i in 1:2)
{
  matplot(t(X[fac==i,]), col=partition[fac==i], type = 'l', main = paste("GOSE",i), xlab="")
}

library(mvtnorm)
library(rgl)
library(car)

clust <- hclust(as.dist(psm), method='ward.D2') # method represents the linkage
View(psm)
View(hclust)

{x11()
plot(clust, main='', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')}


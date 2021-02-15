
remove(list = ls())
#setwd("C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering/Smoothing/Smoothing_def")

source('Cross validation.R')
load('X.RData')

# cut tails
X_1 <- X[,seq(151,1050)]
matplot(t(X_1), type='l')
dim(X_1)
X <- X_1

# set step
step <- 9
n_time    <- dim(X)[2]/step
time.grid <- 1:n_time 
X <- X[, seq(1,dim(X)[2],by = step)]

# spline order
m <- 4

# number of data
n <- dim(X)[1]

### CV for number of basis ----
nbasis <- 6:80


gcv <- matrix(0,n,length(nbasis))  #(gcv)_ij: gcv for datum i smoothed with a j+min_basis number of basis

pb <- txtProgressBar(min = 0, max = length(nbasis), style = 3)

for (i in 1:length(nbasis)){ # Cross validation
  basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=nbasis[i], norder=m)
  gcv[,i] <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)$gcv
  setTxtProgressBar(pb, i)
}

# GCV for each smoothed function vs. number of basis
par(mfrow=c(2,2))
for (j in 1:26) {
  plot(nbasis,gcv[j,], main = j)
  text(nbasis,gcv[j,], labels=30:80)
  optimal <- nbasis[which(min(gcv[j,])==gcv[j,])]
  
}
dev.off()


# mean of GCV over all the functions vs. number of basis
mean_gcv<-colMeans(gcv)
#x11()
par(mfrow=c(1,1))
png(file = paste0("pics/GCV.png"), width = 8000, height = 5000, units = "px", res = 800)
plot(nbasis,mean_gcv, type='l', xlab="Number of basis", ylab="Average GCV", main="GCV", lwd=2, col="deepskyblue4")
grid()
dev.off()
#text(nbasis,mean_gcv, labels=6:80)
nbasis[which.min(mean_gcv)]


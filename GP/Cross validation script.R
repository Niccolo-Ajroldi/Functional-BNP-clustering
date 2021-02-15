
# spline order
m <- 4

# number of data
n <- dim(X)[1]

### CV for number of basis ----
nbasis <- 6:100

gcv <- matrix(0,n,length(nbasis))  #(gcv)_ij: gcv for datum i smoothed with a j+min_basis number of basis

pb <- txtProgressBar(min = 0, max = length(nbasis), style = 3)

for (i in 1:length(nbasis)){ # Cross validation
  basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=nbasis[i], norder=m)
  gcv[,i] <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)$gcv
  setTxtProgressBar(pb, i)
}

# GCV for each smoothed function vs. number of basis
par(mfrow=c(2,2))
for (j in 1:n) {
  plot(nbasis,gcv[j,], main = j)
  text(nbasis,gcv[j,], labels=30:80)
  optimal <- nbasis[which(min(gcv[j,])==gcv[j,])]
  
}
dev.off()


# mean of GCV over all the functions vs. number of basis
mean_gcv<-colMeans(gcv)
x11()
par(mfrow=c(1,1))
plot(nbasis,mean_gcv)
text(nbasis,mean_gcv, labels=6:80)
nbasis[which.min(mean_gcv)]
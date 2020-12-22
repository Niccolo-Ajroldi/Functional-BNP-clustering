library(fda)

remove(list = ls())

setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
load("data_extraction.RData")

#choose one component
my_data<-f.data$ausxSL

#### LS spline regression without penalization

setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
load("data_extraction.RData")

time<-1:1600
X<-my_data$data

#### One observation: -----

## 4th order spline
#create basis
# Set parameters
m <- 4           # spline order 
nbasis <- 9

basis <- create.bspline.basis(rangeval=c(0,1600), nbasis=nbasis, norder=m)
plot(basis)


X_smoothed_f <- smooth.basis(argvals=time, y=X[1,], fdParobj=basis)
X_smoothed <- eval.fd(time, X_smoothed_f$fd)
X_smoothed1 <- eval.fd(time, X_smoothed_f$fd, Lfd = 1)
X_smoothed2 <- eval.fd(time, X_smoothed_f$fd, Lfd = 2)


x11()
par(mfrow=c(1,2))
plot(time, X[1,], type = "l", main = "Original datum", ylim= c(-250,250), ylab='')
plot(time, X_smoothed, type = "l", main = "Smoothed datum", ylim= c(-250,250),ylab='')


x11()
plot(time, X_smoothed, type = "l", main = "Original datum", ylim= c(-250,250), ylab='')
points(time,X[1,],col="green",lwd=2)
#clearly oversmoothing


#CV for number of basis
nbasis <- 6:80
gcv <- numeric(length(nbasis))
pb <- txtProgressBar(min = 0, max = length(nbasis), style = 3)
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(rangeval=c(0,1600), nbasis[i], m)
  gcv[i] <- smooth.basis(time, X[1,], basis)$gcv
  setTxtProgressBar(pb, i)
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
text(nbasis,gcv, labels=6:80)
nbasis[which.min(gcv)] #76
#BUT we observe a "stable" elbow around 45

# 6:60
which(nbasis==60)
par(mfrow=c(1,1))
plot(nbasis[1:55],gcv[1:55])
text(nbasis[1:55],gcv[1:55], labels=6:60)
nbasis[which.min(gcv[1:55])] #56

# 6:55
par(mfrow=c(1,1))
plot(nbasis[1:50],gcv[1:50])
text(nbasis[1:50],gcv[1:50], labels=6:55)
nbasis[which.min(gcv[1:50])] #49


#we perform smoothing with the number of basis = 49
nbasis<-49
basis <- create.bspline.basis(rangeval=c(0,1600), nbasis=nbasis, norder=m)
plot(basis)

X_smoothed_f <- smooth.basis(argvals=time, y=X[1,], fdParobj=basis)
X_smoothed <- eval.fd(time, X_smoothed_f$fd)
X_smoothed1 <- eval.fd(time, X_smoothed_f$fd, Lfd = 1)
X_smoothed2 <- eval.fd(time, X_smoothed_f$fd, Lfd = 2)


x11()
par(mfrow=c(1,2))
plot(time, X[1,], type = "l", main = "Original datum", ylim= c(-250,250), ylab='')
plot(time, X_smoothed, type = "l", main = "Smoothed datum", ylim= c(-250,250),ylab='')


x11()
plot(time, X_smoothed, type = "l", main = "Original datum", ylim= c(-250,250), ylab='')
points(time,X[1,],col="green",lwd=2)


## 3rd order spline --> we do not need the 2nd derivative
#create basis
# Set parameters
m <- 3          # spline order 



#CV for number of basis
nbasis <- 6:80
gcv <- numeric(length(nbasis))
pb <- txtProgressBar(min = 0, max = length(nbasis), style = 3)
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(rangeval=c(0,1600), nbasis[i], m)
  gcv[i] <- smooth.basis(time, X[1,], basis)$gcv
  setTxtProgressBar(pb, i)
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
text(nbasis,gcv, labels=6:80)
nbasis[which.min(gcv)] #79
#BUT STILL "stable" elbow around 45

# 6:60
par(mfrow=c(1,1))
plot(nbasis[1:55],gcv[1:55])
text(nbasis[1:55],gcv[1:55], labels=6:60)
nbasis[which.min(gcv[1:55])] #59

# 6:55
par(mfrow=c(1,1))
plot(nbasis[1:50],gcv[1:50])
text(nbasis[1:50],gcv[1:50], labels=6:55)
nbasis[which.min(gcv[1:50])] #still 49

#smoothing with 52 basis
nbasis<-52
basis <- create.bspline.basis(rangeval=c(0,1600), nbasis=nbasis, norder=m)
plot(basis)

X_smoothed_f <- smooth.basis(argvals=time, y=X[1,], fdParobj=basis)
X_smoothed <- eval.fd(time, X_smoothed_f$fd)
X_smoothed1 <- eval.fd(time, X_smoothed_f$fd, Lfd = 1)
X_smoothed2 <- eval.fd(time, X_smoothed_f$fd, Lfd = 2)


x11()
par(mfrow=c(1,2))
plot(time, X[1,], type = "l", main = "Original datum", ylim= c(-250,250), ylab='')
plot(time, X_smoothed, type = "l", main = "Smoothed datum", ylim= c(-250,250),ylab='')


x11()
plot(time, X_smoothed, type = "l", main = "Original vs smoothed datum", ylim= c(-250,250), ylab='')
points(time,X[1,],col="green",lwd=2)
#similar




#### All the observation ----

## 3rd order spline
# Set parameters
m <- 3          # spline order 



#CV for number of basis
nbasis <- 6:80
gcv <- matrix(0,26,length(nbasis)) 
#(gcv)_ij: gcv for datum i smoothed with a j+6 number of basis
pb <- txtProgressBar(min = 0, max = 26, style = 3)
for(count in 1:26){ #choose a function
  for (i in 1:length(nbasis)){
    basis <- create.bspline.basis(rangeval=c(0,1600), nbasis[i], m)
    gcv[count,i] <- smooth.basis(time, X[count,], basis)$gcv
  }
  setTxtProgressBar(pb, count)
}

mean_gcv<-colMeans(gcv)
par(mfrow=c(1,1))
plot(nbasis,mean_gcv)
text(nbasis,mean_gcv, labels=6:80)
nbasis[which.min(mean_gcv)] #80 BUT elbow around 45


cut_nbasis<-30:80
which(nbasis==30)
length(nbasis)
cut_mean_cvg<-mean_gcv[25:75] 


#for low numbers of basis very high values
# 30:80
par(mfrow=c(1,1))
plot(nbasis[25:75],mean_gcv[25:75])
text(nbasis[25:75],mean_gcv[25:75], labels=30:80)
nbasis[which(min(mean_gcv[25:75])==mean_gcv)] #80


# 40:80
par(mfrow=c(1,1))
plot(nbasis[35:75],mean_gcv[35:75])
text(nbasis[35:75],mean_gcv[35:75], labels=40:80)
nbasis[which(min(mean_gcv[35:75])==mean_gcv)] #80


#try with 60
nbasis<-60
basis <- create.bspline.basis(rangeval=c(0,1600), nbasis=nbasis, norder=m)
plot(basis)

#you can load file 'smooth_60b_nopenalization.RData' and go to plots
X_smoothed <- matrix(0,26,1600) #(X_smoothed)_ij: X_smoothed_i(j), i-th function smoothed at time j
X_smoothed1 <- matrix(0,26,1600) #as above but first derivative
pb <- txtProgressBar(min = 0, max = 26, style = 3)
for (count in 1:26) {
  X_smoothed_f <-smooth.basis(argvals=time, y=X[count,], fdParobj=basis)
  X_smoothed[count,] <- eval.fd(time, X_smoothed_f$fd)
  X_smoothed1[count,] <- eval.fd(time, X_smoothed_f$fd, Lfd = 1)
  setTxtProgressBar(pb, count)
}


#plots
for (count in 1:26) {
  x11()
  par(mfrow=c(1,2))
  plot(time, X[count,], type = "l", main = paste("Original datum",count), ylim= c(-250,250), ylab='')
  plot(time, X_smoothed[count,], type = "l", main = paste("Smoothed datum",count), ylim= c(-250,250),ylab='')
  
}

graphics.off()

for(count in 1:26){
  x11()
  plot(time, X_smoothed[count,], type = "l", main = paste("Original vs smoothed datum",count), ylim= c(-250,250), ylab='')
  points(time,X[count,],col="green",lwd=2)
}



graphics.off()


#first derivatives: we do not need them. Should we look at them for diagnostics? 
#Anyway fanno entrambe letteralmente schifo
rappincX <- (X[,3:1600]-X[,1:(1600-2)])/(time[3:1600]-time[1:(1600-2)])
for (count in 1:26) {
  x11()
  par(mfrow=c(1,2))
  plot(time[2:(1600-1)], rappincX[count,], type = "l", main = paste("Finite differences",count), ylim= c(-250,250),ylab='')
  plot(time, X_smoothed1[count,], type = "l", main = paste("Smoothed datum 1st derivative",count), ylim= c(-250,250),ylab='')
}

graphics.off()



save(X,X_smoothed,X_smoothed1,X_smoothed_f,rappincX, file="smooth_60b_nopenalization.RData")






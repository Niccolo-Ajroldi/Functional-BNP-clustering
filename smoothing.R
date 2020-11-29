library(fda)
library(fda.usc)
library(knitr)
library(fields)
library(compositions)

setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
load("data_extraction.RData")

# extract only one componet per patient

my_data = f.data$ausxSL

x11()
plot(my_data$argvals, my_data$data[1,], type = "l", ylim= c(-250,250), xlab=expression(paste("t","  (",mu,"s)")), ylab= "y(t)", lwd = 2)

for(i in 2:26){
  lines(my_data$argvals, my_data$data[i,], col=i)
}

curves <- my_data$data
t_ax <- as.matrix(1:1600)


# Plot of original functions
x11()
matplot(t_ax,t(curves), type='l', xlab='x', ylab='orig.func')


# First Approach - Spline Regression-----------------------------------------

## Do spline regression (n_knots = 1600) with B-splines of order 4 and penalizing the 1st derivative
## Does it make sense? I know that my ideal curve should have two peaks. Would it be better to penalize
## 2nd derivative in order to have fewer "spurious" peaks?)

# Penalize first derivative

n <- dim(curves)[1]
m <- 4

pen_smooth0 <- matrix(0, dim(curves)[1], dim(curves)[2])
pen_smooth1 <- matrix(0, dim(curves)[1], dim(curves)[2])
pen_smooth2 <- matrix(0, dim(curves)[1], dim(curves)[2])


basis <- create.bspline.basis(breaks = t_ax, norder=m)

functionalPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=1e3)

# WARNING! It takes forever, go down and load "smoothed_data1e3.RData"

# test on only one function
penal <- smooth.basis(t_ax, curves[1,], functionalPar)
penal0 <- eval.fd(t_ax, penal$fd, Lfd = 0)
penal1 <- eval.fd(t_ax, penal$fd, Lfd = 1)
penal2 <- eval.fd(t_ax, penal$fd, Lfd = 2)

x11()
par(mfrow = c(2,2))
plot(t_ax, curves[1,], type = "l", main = "Original Data", ylim= c(-250,250))
plot(t_ax, penal0, type = "l", main = "Smoothed Data", ylim= c(-250,250))
plot(t_ax, penal1, type = "l", main = "Smoothed Data - First Derivative", ylim= c(-250,250))
plot(t_ax, penal2, type = "l", main = "Smoothed Data - Second Derivative", ylim= c(-250,250))

# Smoothing the dataset
for (i in 1:n) {
  penalized = smooth.basis(t_ax, curves[i,], functionalPar)
  p0 = eval.fd(t_ax, penalized$fd, Lfd=0)
  pen_smooth0[i,] = t(p0)
  p1 = eval.fd(t_ax, penalized$fd, Lfd=1)
  pen_smooth1[i,] = t(p1)
  p2 = eval.fd(t_ax, penalized$fd, Lfd=2)
  pen_smooth2[i,] = t(p2)
}


load("smooth_1600b_1D_1e3.RData")

# compare original and smoothed data
x11()
par(mfrow=c(2,2))
matplot(t_ax, t(curves), type = 'l', main = 'Original Data', ylim = c(-250,250))
matplot(t_ax, t(pen_smooth0), type = 'l', main = 'Smoothed Data', ylim = c(-250,250))
matplot(t_ax, t(pen_smooth1), type = 'l', main = 'Original Data - First Derivative')
matplot(t_ax, t(pen_smooth2), type = 'l', main = 'Smoothed Data - Second Derivative')

#for one curve
i<-1
x11()
par(mfrow=c(2,2))
matplot(t_ax, curves[i,], type = 'l', main = 'Original Data', ylim = c(-250,250))
matplot(t_ax, pen_smooth0[i,], type = 'l', main = 'Smoothed Data', ylim = c(-250,250))
matplot(t_ax, pen_smooth1[i,], type = 'l', main = 'Original Data - First Derivative')
matplot(t_ax, pen_smooth2[i,], type = 'l', main = 'Smoothed Data - Second Derivative')

x11()
par(mfrow=c(1,2))
matplot(t_ax, t(curves), type = 'l', main = 'Original Data')
matplot(t_ax, t(pen_smooth0), type = 'l', main = 'Smoothed Data')


# Penalize 2nd derivative, lambda = 1e3

functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e3)

# Smoothing the dataset
for (i in 1:n) {
  penalized = smooth.basis(t_ax, curves[i,], functionalPar)
  p0 = eval.fd(t_ax, penalized$fd, Lfd=0)
  pen_smooth0[i,] = t(p0)
  p1 = eval.fd(t_ax, penalized$fd, Lfd=1)
  pen_smooth1[i,] = t(p1)
  p2 = eval.fd(t_ax, penalized$fd, Lfd=2)
  pen_smooth2[i,] = t(p2)
}

#save(curves, pen_smooth0, pen_smooth1, pen_smooth2, file = "smooth_1600b_2D_1e3.RData")
load("smooth_1600b_2D_1e3.RData")

# compare original and smoothed data
x11()
par(mfrow=c(2,2))
matplot(t_ax, t(curves), type = 'l', main = 'Original Data')
matplot(t_ax, t(pen_smooth0), type = 'l', main = 'Smoothed Data')
matplot(t_ax, t(pen_smooth1), type = 'l', main = 'Original Data - First Derivative')
matplot(t_ax, t(pen_smooth2), type = 'l', main = 'Smoothed Data - Second Derivative')

#for one curve
i<-1
x11()
par(mfrow=c(2,2))
matplot(t_ax, curves[i,], type = 'l', main = 'Original Data', ylim = c(-250,250))
matplot(t_ax, pen_smooth0[i,], type = 'l', main = 'Smoothed Data', ylim = c(-250,250))
matplot(t_ax, pen_smooth1[i,], type = 'l', main = 'Original Data - First Derivative')
matplot(t_ax, pen_smooth2[i,], type = 'l', main = 'Smoothed Data - Second Derivative')


x11()
par(mfrow=c(1,2))
matplot(t_ax, t(curves), type = 'l', main = 'Original Data')
matplot(t_ax, t(pen_smooth0), type = 'l', main = 'Smoothed Data')

#Boh mi sembra uguale a quando penalizzo la derivata prima



# Smoothing with 500 knots ------------------------------------------------

# I don't use the points as knots. Instead, I decrease the number of knots to 500
# smooth.basis command becomes faster

range <- range(t_ax)
nbasis_ <- 500
basis <- create.bspline.basis(range, nbasis = nbasis_, norder=m)

functionalPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=1e3)

# Smoothing the dataset --> jump down to load .RData
for (i in 1:n) {
  penalized = smooth.basis(t_ax, curves[i,], functionalPar)
  p0 = eval.fd(t_ax, penalized$fd, Lfd=0)
  pen_smooth0[i,] = t(p0)
  p1 = eval.fd(t_ax, penalized$fd, Lfd=1)
  pen_smooth1[i,] = t(p1)
  p2 = eval.fd(t_ax, penalized$fd, Lfd=2)
  pen_smooth2[i,] = t(p2)
}

#save(curves, pen_smooth0, pen_smooth1, pen_smooth2, file = "smooth_500b_1D_1e3.RData")
load("smooth_500b_1D_1e3.RData")

x11()
par(mfrow=c(2,2))
matplot(t_ax, t(curves), type = 'l', main = 'Original Data')
matplot(t_ax, t(pen_smooth0), type = 'l', main = 'Smoothed Data')
matplot(t_ax, t(pen_smooth1), type = 'l', main = 'Smoothed Data - First Derivative')
matplot(t_ax, t(pen_smooth2), type = 'l', main = 'Smoothed Data - Second Derivative')


#for one curve
i<-1
x11()
par(mfrow=c(2,2))
matplot(t_ax, curves[i,], type = 'l', main = 'Original Data', ylim = c(-250,250))
matplot(t_ax, pen_smooth0[i,], type = 'l', main = 'Smoothed Data', ylim = c(-250,250))
matplot(t_ax, pen_smooth1[i,], type = 'l', main = 'Original Data - First Derivative')
matplot(t_ax, pen_smooth2[i,], type = 'l', main = 'Smoothed Data - Second Derivative')


# Penalize 2nd derivative

basis <- create.bspline.basis(range, nbasis = nbasis_, norder=m)

functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e3)

# Smoothing the dataset
for (i in 1:n) {
  penalized = smooth.basis(t_ax, curves[i,], functionalPar)
  p0 = eval.fd(t_ax, penalized$fd, Lfd=0)
  pen_smooth0[i,] = t(p0)
  p1 = eval.fd(t_ax, penalized$fd, Lfd=1)
  pen_smooth1[i,] = t(p1)
  p2 = eval.fd(t_ax, penalized$fd, Lfd=2)
  pen_smooth2[i,] = t(p2)
}

#save(curves, pen_smooth0, pen_smooth1, pen_smooth2, file = "smooth_500b_2D_1e3.RData")
load("smooth_500b_2D_1e3.RData")

x11()
par(mfrow=c(2,2))
matplot(t_ax, t(curves), type = 'l', main = 'Original Data')
matplot(t_ax, t(pen_smooth0), type = 'l', main = 'Smoothed Data')
matplot(t_ax, t(pen_smooth1), type = 'l', main = 'Smoothed Data - First Derivative')
matplot(t_ax, t(pen_smooth2), type = 'l', main = 'Smoothed Data - Second Derivative') 


#for one curve
i<-1
x11()
par(mfrow=c(2,2))
matplot(t_ax, curves[i,], type = 'l', main = 'Original Data', ylim = c(-250,250))
matplot(t_ax, pen_smooth0[i,], type = 'l', main = 'Smoothed Data', ylim = c(-250,250))
matplot(t_ax, pen_smooth1[i,], type = 'l', main = 'Original Data - First Derivative')
matplot(t_ax, pen_smooth2[i,], type = 'l', main = 'Smoothed Data - Second Derivative')




# (FAILED) temptative of CV for the choice of lambda ----------------------


# Cross validation for the choice of lambda
# 
# I do it on a single curve, the one with the steepest peak
curve <- curves[which(curves[,500] > 700),]

nbasis_ <- 300
m <- 4 

basis = create.bspline.basis(range, nbasis = nbasis_, norder=m)

lambda=c(1,1e1,1e2,1e3,1e4,1e5)
gcv=numeric(length(lambda))
for (i in 1:length(lambda))
{
  functionalPar = fdPar(fdobj=basis, Lfdobj=1, lambda=lambda[i])  
  gcv[i] = smooth.basis(t_ax, curve, functionalPar)$gcv
}



plot(log10(lambda),gcv)
lambda[which.min(gcv)]

# gcv mi dice che dovrei scegliere lambda piccolo, ma così non faccio undersmoothing?

functionalPar = fdPar(fdobj=basis, Lfdobj=1, lambda=1e2)
penalized <- smooth.basis(t_ax, curve, functionalPar)
smoothed_curve = eval.fd(t_ax, penalized$fd, Lfd=0)

x11()
par(mfrow = c(1,2))
plot(t_ax, curve, type="l", col = "red", main = "Original Curve")
plot(t_ax, smoothed_curve, type="l", col = "blue", main = "Smoothe Curve")




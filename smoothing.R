library(fda)
library(fda.usc)
library(knitr)
library(fields)
library(compositions)

setwd("C:/Users/Teresa Bortolotti/Desktop/ProgettoBayesiana")
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

#There seem to be three clusters of functions: in fact, one of them is obtained by warping the abscissas.


# Plot of original functions
x11()
matplot(t_ax,t(curves), type='l', xlab='x', ylab='orig.func')

n <- dim(curves)[1]
m <- 4

pen_smooth0 <- matrix(0, dim(curves)[1], dim(curves)[2])
pen_smooth1 <- matrix(0, dim(curves)[1], dim(curves)[2])
pen_smooth2 <- matrix(0, dim(curves)[1], dim(curves)[2])


basis <- create.bspline.basis(t_ax, norder=m)

functionalPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=1e3)

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


# compare original and smoothed data
x11()
par(mfrow=c(1,2))
matplot(t_ax, t(curves), type = 'l', main = 'Original Data', ylim = c(-250,250))
matplot(t_ax, t(pen_smooth0), type = 'l', main = 'Smoothed Data', ylim = c(-250,250))


x11()
par(mfrow=c(1,2))
matplot(t_ax, t(curves), type = 'l', main = 'Original Data')
matplot(t_ax, t(pen_smooth0), type = 'l', main = 'Smoothed Data')


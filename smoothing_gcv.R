# data ----------------------

setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
load("data_extraction.RData")

# extract only one componet per patient

my_data = f.data$ausxSL
curves <- my_data$data
t_ax <- as.matrix(1:1600)
n <- dim(curves)[1]

# Cross validation for the choice of lambda ----------------------

nbasis <- 50
m <- 4

# generate a basis
basis <- create.bspline.basis(rangeval=range(t_ax),
                             nbasis = nbasis, 
                             norder=m)

# grid of lambda
lambda <- 10^seq(-10,10)
gcv <- matrix(0,length(lambda),n)

for (i in 1:length(lambda))
{
  print(i)
  functionalPar = fdPar(fdobj=basis, Lfdobj=1, lambda=lambda[i])  
  gcv[i,] = smooth.basis(t_ax, t(curves), functionalPar)$gcv # see smooth.basis->argument y per il trasposto
}

# plot GCV for all the curves
x11()
matplot(log10(lambda), gcv, type = 'l', lwd=2, xlab="log10(lambda)", ylab="GCV error", main="GCV(lambda) error for different curves")

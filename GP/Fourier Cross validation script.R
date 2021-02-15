
# number of data
n <- dim(X)[1]

# Fourier basis
basis <- create.fourier.basis(rangeval=range(time.grid), nbasis=7)

# coefficients
beta <- Data2fd(y=t(X), argvals=time.grid, basisobj=basis)$coefs

# plot smoothed data
basis.t <- t(eval.basis(time.grid, basis))
X_smooth <- t(beta) %*% basis.t
x11()
matplot(time.grid, t(X_smooth), type='l', col=col)



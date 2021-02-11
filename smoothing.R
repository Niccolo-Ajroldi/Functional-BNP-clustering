
smoothing <- function(X, step=1, nbasis, spline_order=4){
  
  n_time    <- dim(X)[2]/step
  time.grid <- 1:n_time 
  #time.grid <- seq(1,dim(X)[2],by=step) 
  
  # use steps
  X <- X[, seq(1,dim(X)[2],by=step)] # matrix n x time.grid, HO TOLTO ANCHE LA 24

  # basis 
  basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=nbasis, norder=spline_order)
  
  # smoothing
  X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)
  
  # save coefficients
  beta <- t(X_smoothed_f$fd$coefs)
  
  # plot smoothed data
  basis.t <- t(eval.basis(time.grid, basis))
  X_smooth <- beta %*% basis.t
  matplot(time.grid, t(X_smooth), type='l')
  
  smoothing_parameters <- list('step' = step,
                               'number_basis' = nbasis,
                               'spline_order' = spline_order)
  
  return(list('basis' = basis,
              'beta'= beta,
              'time.grid' = time.grid,
              'X' = X,
              'smoothing_parameters' = smoothing_parameters))
}


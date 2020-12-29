
smoothing <- function(X, step, nbasis, spline_order){
  
  n_time    <- 1600/step
  time.grid <- 1:n_time 
  
  # use steps
  X <- X[, seq(1,1600,by=step)] # matrix n x time.grid, HO TOLTO ANCHE LA 24

  # basis 
  basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=nbasis, norder=spline_order)
  
  # smoothing
  X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)
  
  # save coefficients
  beta <- t(X_smoothed_f$fd$coefs)
  
  smoothing_parameters <- list('step' = step,
                               'number_basis' = nbasis,
                               'spline_order' = spline_order)
  
  return(list('basis' = basis,
              'beta'= beta,
              'time.grid' = time.grid,
              'X' = X,
              'smoothing_parameters' = smoothing_parameters))
}
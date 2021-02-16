
#'
#' Smoothing
#' 
#' Given the observed data, the time step, the desired number of basis
#' and the desired spline order, the function constructs the smoothed data,
#' returning a list containing the original data, the parameters of the smoothing,
#' the basis functions and the basis coefficients.
#' 
#' @param X: matrix of observed data (X)_ij = X_i(j), dimensions: (n)x(n_time)
#' @param step: time step that allows us to build a less dense time grid
#'              on which to build the smoothed data
#' @param nbasis: desired number of basis functions
#' @param spline_order: desired order of the splines used to build the basis
#'
#' @return a list containing the following objects:
#'         input values: X, step, nbasis, spline_order
#'         time grid:
#'         basis: functional data object (fda package)
#'         beta: estimated coefficients wrt the basis functions
#'

smoothing <- function(X,
                      step=1,
                      nbasis,
                      spline_order=4)
{
  
  n_time    <- dim(X)[2]/step
  time.grid <- 1:n_time
  
  # use steps
  X <- X[, seq(1,dim(X)[2],by=step)]

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
  
  
  return(list('X' = X,
              'step' = step,
              'time.grid' = time.grid,
              'number_basis' = nbasis,
              'spline_order' = spline_order,
              'basis' = basis,
              'beta'= beta))
}


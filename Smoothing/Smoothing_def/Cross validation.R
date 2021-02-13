###################################################
#### LS spline regression without penalization ####
###################################################


cross_validation_nbasis <- function(X, min_basis = 10, max_basis = 250, step = 1, spline_order = 4){
  

  n_time    <- dim(X)[2]/step
  time.grid <- 1:n_time 

  m <- spline_order
  
  n <- dim(X)[1]
  
  ### CV for number of basis ----
  m <- spline_order
  nbasis <- min_basis:max_basis
  
  
  gcv <- matrix(0,26,length(nbasis))  #(gcv)_ij: gcv for datum i smoothed with a j+ number of basis
  
  
  for (i in 1:length(nbasis)){ # Cross validation
    basis <- create.bspline.basis(rangeval=range(time.grid), nbasis=nbasis, norder=spline_order)
    gcv[,i] <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)
  }
  
}





# input
#  var_phi: desired variance for phi_t 
#  X_obs: matrix (X_obs)_ij = X_i(j)
#  beta: matrix (beta)_ij projection coefficients of datum i on the basis j

new_hyperparameters <- function(var_phi, X_obs, beta, time.grid, scale){
  
  #### phi_t : c,d ####
  ## phi_t ~ IG(c,d) forall t
  mean_phi <- apply(X_obs,MARGIN=2,FUN=var)*scale
  C <- mean_phi^2/var_phi + 2
  D <- mean_phi^3/var_phi + mean_phi
  
  #### mu : m0,Lambda0 ####
  ## mu ~ N_L(m0,Lambda0)
  m0 <- colMeans(beta)
  Lambda0 <- diag(diag(cov(beta)))
  
  desired_prior_values <- c('var_phi' = var_phi)
  
  return (list('C' = C,
               'D' = D,
               'm0' = m0,
               'Lambda0' = Lambda0,
               'desired_prior_values' = desired_prior_values))
}

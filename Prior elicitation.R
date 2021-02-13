
# input
#  var_phi: desired variance for phi_t 
#  X_obs: matrix (X_obs)_ij = X_i(j)
#  beta: matrix (beta)_ij projection coefficients of datum i on the basis j

hyperparameters <- function(mean_phi,var_phi, X_obs, beta, scale=1){
  
  #### phi_t : c,d ####
  ## phi_t ~ IG(c,d) forall t
  #mean_phi <- 100
  #mean_phi <- mean(apply(X_obs,MARGIN=2,FUN=var)) * scale
  
  c <- mean_phi^2/var_phi + 2
  d <- mean_phi^3/var_phi + mean_phi
  if(c<=2 | d<0)
    stop('Parameters c and d not compatible')
  
  #### mu : m0,Lambda0 ####
  ## mu ~ N_L(m0,Lambda0)
  m0 <- colMeans(beta)
  Lambda0 <- cov(beta)
  
  desired_prior_values <- c('var_phi' = var_phi)
  
  return (list('c' = c,
               'd' = d,
               'm0' = m0,
               'Lambda0' = Lambda0,
               'desired_prior_values' = desired_prior_values))
}

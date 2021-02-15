#'
#' Hyperparameters initialization
#' 
#' 
#' Given observed data, basis expansion coefficients and a specified choice of mean_phi and var_phi, 
#' the function initializes a list containing values for the hyperparameters,
#' necessary for the run of FBNP_hyper_alltime sampling.
#' 
#' @param X_obs: matrix of observed data (X_obs)_ij = X_i(j), dimensions: (n)x(n_time)
#' @param beta: matrix (beta)_ij projection coefficients of datum i on the basis j
#' @param mean_phi: desired mean for phi_t 
#' @param var_phi: desired variance for phi_t 
#' 
#' @return a list with initialized hyperparameters for running FBNP of FBNP_hyper.
#' 
hyperparameters_alltime <- function(X_obs,
                                    beta,
                                    var_phi){
  
  #### phi_t : c,d ####
  ## phi_t ~ IG(c,d) forall t
  mean_phi <- apply(X_obs,MARGIN=2,FUN=var)
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
               'var_phi' = var_phi,
               'mean_phi' = mean_phi))
  
}

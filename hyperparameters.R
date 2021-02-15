#'
#' Hyperparameters initialization
#' 
#' 
#' Given observed data, basis expansion coefficients and a specified choice of mean_phi and var_phi, 
#' the function initializes a list containing values for the hyperparameters,
#' necessary for the run of FBNP or FBNP_hyper sampling.
#' 
#' @param X_obs: matrix of observed data (X_obs)_ij = X_i(j), dimensions: (n)x(n_time)
#' @param beta: matrix (beta)_ij projection coefficients of datum i on the basis j
#' @param mean_phi: desired mean for phi_t 
#' @param var_phi: desired variance for phi_t 
#' 
#' @return a list with initialized hyperparameters for running FBNP of FBNP_hyper.
#'         

hyperparameters <- function(X_obs, 
                            beta, 
                            mean_phi,
                            var_phi)
{
  
  #### phi_t : c,d ####
  ## phi_t ~ IG(c,d) forall t
  c <- mean_phi^2/var_phi + 2
  d <- mean_phi^3/var_phi + mean_phi
  if(c<=2 | d<0)
    stop('Parameters c and d not compatible')
  
  #### mu : m0,Lambda0 ####
  ## mu ~ N_L(m0,Lambda0)
  m0 <- colMeans(beta)
  Lambda0 <- cov(beta)
  
  return (list('c' = c,
               'd' = d,
               'm0' = m0,
               'Lambda0' = Lambda0,
               'var_phi' = var_phi,
               'mean_phi' = mean_phi))
}

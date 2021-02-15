
#'
#' PSM
#' 
#' Given the matrix K that contains, for each iteration, the latent partition of the
#' observations, the function returns the Posterior Similarity Matrix
#' 
#' 
#' @param K: matrix of the latent partition (K)_ij: cluster assignment of observation j at iteration i
#'                        
#' @return the Posterior Similarity Matrix of the partition
#'


PSM <- function(K)
{
  n <- ncol(K)
  psm <- matrix(0, ncol = n, nrow = n)
  for(i in 1:n){
    
    for(j in 1:i){
      
      psm[i,j] <- psm[j,i] <- mean(K[,i] == K[,j])
      
    }
  }
  
  return(psm)
  
}
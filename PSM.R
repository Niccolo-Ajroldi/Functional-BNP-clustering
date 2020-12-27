# PSM - Posterior Similarity Matrix

PSM <- function(cls){
  n <- ncol(cls)
  mat_out <- matrix(0, ncol = n, nrow = n)
  for(i in 1:n){
    for(j in 1:i){
      mat_out[i,j] <- mat_out[j,i] <- mean(cls[,i] == cls[,j])
    }
  }
  return(mat_out)
}

#### NUMBER OF CLUSTERS ####

unique_clusters <- apply(K, 1, function(x) length(unique(x)))
coda_import <- as.mcmc(unique_clusters)
summary(coda_import)
geweke.diag(coda_import)
acfplot(coda_import)

#### Posterior Similarity Matrix ####

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

heatmap(PSM(K), Rowv = NA, Colv = NA)

est_part_BINDER <- function(clust, PSM){
  res_loss <- c()
  for(i in 1:nrow(clust)){
    binary_part <- sapply(clust[i,], function(x) as.numeric(x == clust[i,]))
    res_loss[i] <- sum(abs(binary_part - PSM))
  }
  return(clust[which.min(res_loss),])
}


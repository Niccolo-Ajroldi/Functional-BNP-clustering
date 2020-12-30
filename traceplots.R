
#library(installr)
#install.ImageMagick()
library(magick)
library(animation)

names(out)

# output of the algorithm
K            <- out$K
mu_coef_out  <- out$mu_coef_out
sigma2_out   <- out$sigma2_out
probs_j_out  <- out$probs_j_out
probs_ij_out <- out$probs_ij_out

n_iter <- run_parameters$algorithm_parameters[1]
burnin <- run_parameters$algorithm_parameters[2]
nn <- n_iter-burnin

#### Cluster assignment ####------------------------------

x11()
par(mfrow=c(5,5))
for(i in 1:22)
{
  traceplot(as.mcmc(K[,i]), main=paste0("Observation ",i))
}


#### MU ####------------------------------

is(mu_coef_out)
is(mu_coef_out[[1]])
dim(mu_coef_out[[1]])
length(mu_coef_out)

# for now fix a cluster
j = 1

# I would like to obtain a vector for each iteration,
# the vector of values of mu for each iteration

# first I need to retrieve the coefficients in basis for each iter
mu_coef_j <- vector("list", length = n_iter-burnin)
mu_j      <- vector("list", length = n_iter-burnin)
time.grid <- 1:(1600/5)
basis <- smoothing_list$basis
basis.t <- t(eval.basis(time.grid, basis))
for(iter in 1:nn)
{
  # coefficients of basis projection
  mu_coef_j[[iter]] <- mu_coef_out[[iter]][j,]
  # mu evaluated on the time grid
  mu_j[[iter]] <- mu_coef_j[[iter]] %*% basis.t
}

#x11()
#par(mfrow=c(5,5))
#for(i in 1: 25)
#{
#  matplot(t(mu_j[[i]]), type='l')
#}

# traceplot as a gif
saveGIF(
  expr = {
    for(iter in 1:nn){
      matplot(t(mu_j[[iter]]), type='l', main=paste0("Traceplot of mu(t) for cluster 1, iter=",iter) )
    }
  },
  movie.name = "cluster_1_mu.gif",
  interval = 0.2
)

# TODO: fissare gli assi
# magari far vedere solo i clusters finali significativi?
# aggiungere iterazionis



# wraphttps://stackoverflow.com/questions/28810305/r-invalidargument-delay-with-animation-and-ggplot








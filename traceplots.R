
#library(installr)
#install.ImageMagick()
library(magick)
library(animation)
library(coda)
library(pbmcapply)


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
#mu_coef_j <- vector("list", length = n_iter-burnin)
#mu_j      <- vector("list", length = n_iter-burnin)
#time.grid <- 1:(1600/5)
#basis <- smoothing_list$basis
#basis.t <- t(eval.basis(time.grid, basis))
#for(iter in 1:nn)
#{
#  # coefficients of basis projection
#  mu_coef_j[[iter]] <- mu_coef_out[[iter]][j,]
#  # mu evaluated on the time grid
#  mu_j[[iter]] <- mu_coef_j[[iter]] %*% basis.t
#}

########### PLOT 1 KERNELS ###############################################

saveGIF(
  expr = {
    for(iter in 1:nn){
      par(mfrow=c(3,5))
      par(oma=c(0,0,2,0))
      j=1
      # coefficients of basis projection
      mu_coef_j <- mu_coef_out[[iter]][j,]
      # mu evaluated on the time grid
      mu_j <- mu_coef_j %*% basis.t
      # plot
      matplot(t(mu_j), type='l', main=paste0("Kernel ",j), xlab="Iterations", ylim=c(-300,300))
      #mtext(paste0("Traceplots of mu(t), iteration = ", iter))
      title(paste0("mu(t), iteration ", iter), outer = TRUE)
    }
  },
  movie.name = "cluster_1_15_mu.gif",
  interval = 0.2, 
  ani.width = 1300,
  ani.height = 700
)


########### PLOT 15 KERNELS ###############################################

saveGIF(
  expr = {
    for(iter in 1:nn){
      layout(1)
      par(oma=c(0,0,2,0))
      for(j in 1:15)
      {
        # coefficients of basis projection
        mu_coef_j <- mu_coef_out[[iter]][j,]
        # mu evaluated on the time grid
        mu_j <- mu_coef_j %*% basis.t
        # plot
        matplot(t(mu_j), type='l', main=paste0("Kernel ",j), xlab="Iterations", ylim=c(-300,300))
      }
      #mtext(paste0("Traceplots of mu(t), iteration = ", iter))
      title(paste0("mu(t), iteration ", iter), outer = TRUE)
    }
  },
  movie.name = "cluster_1_15_mu.gif",
  interval = 0.2, 
  ani.width = 1300,
  ani.height = 700
)

############### PLOT ALL KERNELS 15 each time ########################################

setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

# save old directory to get back there at the end
dir.current <- getwd()

# if there is not a Traceplot directory throw an error and stop
if(!dir.exists("Traceplots"))
  stop("Traceplots directory does not exists")

# name of directory where I will put plots, I use current time in the name
new.dir <- paste0("Traceplots",format(Sys.time(),'_%Y%m%d_%H%M%S'))
new.dir <- paste0(dir.current,"/","Traceplots/",new.dir)

# create such directory and go there
dir.create(new.dir)
setwd(new.dir) 


# create 10 folders, named with the kernels in each folder
for(yo in 1:10)
{
  # kernels in this folder
  kernelz.to.plot <- (1+15*(yo-1)):(15*yo)
  
  # minimum and maximum kernel present here
  kernel.lower <- min(kernelz.to.plot)
  kernel.upper <- max(kernelz.to.plot)
  
  # create folder
  dir.create(paste0("Kernel ", kernel.lower, "-", kernel.upper))
}



#######################################################################################################
## MU ##

time.grid <- 1:(1600/5) # TODO: step hereee ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
basis <- smoothing_list$basis
basis.t <- t(eval.basis(time.grid, basis))


# for the 10 "blocks"
for(yo in 1:10)
{
  # I plot 15 kernels in each gif
  kernelz.to.plot <- (1+15*(yo-1)):(15*yo)

  # minimum and maximum kernel present here
  kernel.lower <- min(kernelz.to.plot)
  kernel.upper <- max(kernelz.to.plot)
  
  # go to the folder
  setwd(new.dir)
  setwd(paste0("Kernel ", kernel.lower, "-", kernel.upper))
  
  # traceplot as a gif
  saveGIF(
    expr = {
      for(iter in 1:nn){
        par(mfrow=c(3,5))
        par(oma=c(0,0,2,0))
        for(j in kernelz.to.plot)
        {
          # coefficients of basis projection
          mu_coef_j <- mu_coef_out[[iter]][j,]
          # mu evaluated on the time grid
          mu_j <- mu_coef_j %*% basis.t
          # plot
          matplot(t(mu_j), type='l', main=paste0("Kernel ",j), xlab="Iterations", ylim=c(-300,300))
        }
        #mtext(paste0("Traceplots of mu(t), iteration = ", iter))
        title(paste0("mu(t), iteration ", iter), outer = TRUE)
      }
    },
    movie.name = paste0("Kernel_",kernel.lower,"_to_",kernel.upper,"_mu.gif"),
    interval = 0.2, 
    ani.width = 1300,
    ani.height = 700
  )
  
  # go back to current traceplot folder
  setwd(new.dir)
  
  # ProgressBar
  setTxtProgressBar(pb, yo)
}


# TODO: fissare gli assi
# magari far vedere solo i clusters finali significativi?
# aggiungere iterazionis

# wraphttps://stackoverflow.com/questions/28810305/r-invalidargument-delay-with-animation-and-ggplot


#### SIGMA ################################################################

# matrix, Mx(n_iter-burniin), row j contains the explored values for sigma2 of kernel j
sigma2.out.matrix <- matrix(unlist(sigma2_out), nrow=150, ncol=nn, byrow=FALSE)
#max(sigma2out.matrix[,16] -  sigma2_out[[16]])

# go back to current traceplot folder
setwd(new.dir)

x11(width=1300, height=700)
par(mfrow=c(3,5))
par(oma=c(0,0,2,0))
for(j in 1:15)
{
  traceplot(as.mcmc(sigma2.out.matrix[j,]), main=paste0("Kernel ",j) )#, ylim=c(1e-17,1e-18))
  abline(h=0, lty=2, col="red")
}
title(paste0("sigma^2"), outer = TRUE)

#### PHI ################################################################

# qui per ogni cluster devo fare 1600 traceplots!!!

setwd(dir.current)

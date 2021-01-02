
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

#library(installr)
#install.ImageMagick()
library(coda)
library(pbmcapply)

fun_traceplots <- function (out,
                            smoothing_list,
                            run_parameters,
                            blocks=1,
                            time=0.1,
                            mu_gif=FALSE,
                            iter_step=10) # gif con 100 ci mette un fracco (10 min)
                                          # modificare per vedere solo 10 cluster per volta
  
{
  
  if(mu_gif)
  {
    library(magick)
    library(animation)
  }
    
  basis <- smoothing_list$basis
  time.grid <- smoothing_list$time.grid
  #step <- smoothing_list$smoothing_parameters$step
  
  K            <- out$K
  mu_coef_out  <- out$mu_coef_out
  sigma2_out   <- out$sigma2_out
  probs_j_out  <- out$probs_j_out
  probs_ij_out <- out$probs_ij_out
  
  n_iter <- run_parameters$algorithm_parameters[1]
  burnin <- run_parameters$algorithm_parameters[2]
  thin   <- run_parameters$algorithm_parameters[3]
  M      <- run_parameters$algorithm_parameters[4]
  mass   <- run_parameters$algorithm_parameters[5]
  
  nn <- n_iter-burnin
  
  #### Traceplot Directory ####------------------------------------------------------------
  
  # save old directory to get back there at the end
  dir.current <- getwd()
  
  # if there is not a Result directory throw an error and stop
  if(!dir.exists("Results"))
    stop("Results directory does not exists")
  
  # if there is not a Traceplot directory throws an error and stop
  if(!dir.exists("Results/Traceplots"))
    stop("Results/Traceplots directory does not exists")
  
  # name of directory where I will put plots, I use current time in the name
  new.dir <- paste0("Traceplots",format(Sys.time(),'_%Y%m%d_%H%M%S'))
  new.dir <- paste0(dir.current,"/Results/Traceplots/",new.dir)
  
  # create such directory and go there
  dir.create(new.dir)
  setwd(new.dir) 
  
  # save RData
  save(out, run_parameters, file = "Output.RData")
  
  cat(
    "ALGORITHM\n",
    "Iter:   ", n_iter, "\n",
    "Burnin: ", burnin, "\n",
    #"Thinning:  ", thin, "\n",
    "M:      ", M, "\n",
    "Mass:   ", mass, "\n\n",
    
    "PRIOR\n",
    "var_sigma: ", run_parameters$prior_parameters$desired_prior_values[[1]], "\n",
    "var_phi:   ", run_parameters$prior_parameters$desired_prior_values[[2]], "\n",
    "a:         ", run_parameters$prior_parameters$a, "\n",
    "b:         ", run_parameters$prior_parameters$b, "\n",
    "c:         ", run_parameters$prior_parameters$c, "\n",
    "d:         ", run_parameters$prior_parameters$d, "\n\n",
    
    "SMOOTHING\n",
    "step:       ", run_parameters$smoothing_parameters$step, "\n",
    "L:          ", run_parameters$smoothing_parameters$number_basis, "\n",
    "rescale:    ", run_parameters$smoothing_parameters$rescale_parameter, "\n",
    "eliminated: ", run_parameters$smoothing_parameters$observation_eliminated,
    
    file="Output.txt"
  )
  
  #### Cluster assignment ####-----------------------------------------------------------
  
  X11(width=1300, height=700)
  par(mfrow=c(5,5))
  par(oma=c(0,0,2,0))
  for(i in 1:22)
  {
    traceplot(as.mcmc(K[,i]), main=paste0("Observation ",i), ylim=c(0,150))
  }
  title("Cluster assignment", outer = TRUE)
  savePlot(filename = "Cluster_traceplot",
           type = c("png"),
           device = dev.cur())
  
  
  #### Sigma2 ####------------------------------------------------------------------------
  
  # matrix, Mx(n_iter-burniin), row j contains the explored values for sigma2 of kernel j
  sigma2.out.matrix <- matrix(unlist(sigma2_out), nrow=150, ncol=nn, byrow=FALSE)
  #max(sigma2out.matrix[,16] -  sigma2_out[[16]])
  
  # for "blocks"
  for(yo in blocks)
  {
    # I plot 10 kernels in each gif
    kernelz.to.plot <- (1+10*(yo-1)):(10*yo)
    
    # minimum and maximum kernel present here
    kernel.lower <- min(kernelz.to.plot)
    kernel.upper <- max(kernelz.to.plot)
    
    x11(width=1300, height=700)
    par(mfrow=c(2,5))
    par(oma=c(0,0,2,0))
    for(j in kernelz.to.plot)
    {
      traceplot(as.mcmc(sigma2.out.matrix[j,]), main=paste0("Kernel ",j) )#, ylim=c(1e-17,1e-18))
      abline(h=0, lty=2, col="red")
    }
    title(paste0("sigma^2"), outer = TRUE)
    savePlot(filename = paste0("sigma_",kernel.lower,"_",kernel.upper),
             type = c("png"),
             device = dev.cur())
  }
  
  
  #### MU ####----------------------------------------------------------------------------------
  
  if(mu_gif)
  {
    #is(mu_coef_out)
    #is(mu_coef_out[[1]])
    #dim(mu_coef_out[[1]])
    #length(mu_coef_out)
    
    basis <- smoothing_list$basis
    basis.t <- t(eval.basis(time.grid, basis))
    
    # for "blocks"
    for(yo in blocks)
    {
      # I plot 10 kernels in each gif
      kernelz.to.plot <- (1+10*(yo-1)):(10*yo)
      
      # minimum and maximum kernel present here
      kernel.lower <- min(kernelz.to.plot)
      kernel.upper <- max(kernelz.to.plot)
      
      print(paste0("Making a gif for mu(t), kernels: ",kernel.lower,"_",kernel.upper))
      pb <- progressBar(0, max = nn, initial = 0, style = "ETA")
      
      # traceplot as a gif
      saveGIF(
        expr = {
          for(iter in seq(1,nn,by=iter_step))
          {
            par(mfrow=c(2,5))
            par(oma=c(0,0,2,0))
            for(j in kernelz.to.plot)
            {
              # coefficients of basis projection
              mu_coef_j <- mu_coef_out[[iter]][j,]
              # mu evaluated on the time grid
              mu_j <- mu_coef_j %*% basis.t
              # plot
              matplot(t(mu_j), type='l', main=paste0("Kernel ",j), xlab="Iterations", ylim=c(-100,100))
            }
            #mtext(paste0("Traceplots of mu(t), iteration = ", iter))
            title(paste0("mu(t), iteration ", iter), outer = TRUE)
            
            setTxtProgressBar(pb, iter)
          }
        },
        movie.name = paste0("mu_",kernel.lower,"_",kernel.upper,".gif"),
        interval = 0.1, 
        ani.width = 1300,
        ani.height = 700
      )
      
    }
    
  }
  
  # go back to current traceplot folder
  setwd(dir.current)
  
  
}
  
  
  
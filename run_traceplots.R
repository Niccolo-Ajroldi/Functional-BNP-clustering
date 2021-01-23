
# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

#load("Results/Nico_M50_4_01.RData")

source("fun_traceplots.R")

# save old directory to get back there at the end
dir.current <- getwd()

# call
fun_traceplots (out,
                smoothing_list,
                run_parameters,
                blocks=1,    # each "block" contains 10 kernels, change this in order to see a different gif for mu(t)
                time=0.15,   # time between images in the gif
                mu_gif=FALSE,
                falzo=FALSE  # change this if you are using falzo (different scale of mu(t))
                )

# go to previous directory
setwd(dir.current)


graphics.off()



##-----------------------------------------------
K            <- out$K
mu_coef_out  <- out$mu_coef_out
sigma2_out   <- out$sigma2_out
probs_j_out  <- out$probs_j_out
probs_ij_out <- out$probs_ij_out
phi_out      <- out$phi_out

n         <- dim(smoothing_list$X)[1]
basis     <- smoothing_list$basis
time.grid <- smoothing_list$time.grid
n_time    <- length(time.grid)

n_iter <- run_parameters$algorithm_parameters[[1]]
burnin <- run_parameters$algorithm_parameters[[2]]
thin   <- run_parameters$algorithm_parameters[[3]]
M      <- run_parameters$algorithm_parameters[[4]]
mass   <- run_parameters$algorithm_parameters[[5]]

# number of iterations after burnin
nn <- n_iter-burnin

# show a gif once every iter_step
iter_step = nn/50

phi <- phi_out


for(i in 1:nn)
{
  
  
}

length(phi) # nn
dim(phi[[1]]) # M x n_time

phi_unlist <- unlist(phi)
length(phi_unlist)
M*n_time*nn

cluster <- 1
time <- 1

phi_cluster1_time1 <- phi_unlist[seq(1,length(phi_unlist),by=M*n_time)]
length(phi_cluster1_time1) # torna brooo
traceplot( as.mcmc(phi_cluster1_time1), main=paste0("Kernel ",1,", t=1") )

phi_cluster1_time2 <- phi_unlist[seq(2,length(phi_unlist),by=M*n_time)]
length(phi_cluster1_time1) # torna brooo
traceplot( as.mcmc(phi_cluster1_time2), main=paste0("Kernel ",1,", t=1") )

cluster <- 1
time <- 2
  
phi_cluster_j_time_g <- phi_unlist[seq(time + (cluster-1)*M , length(phi_unlist), by=M*n_time)]
length(phi_cluster1_time1) # torna brooo
traceplot( as.mcmc(phi_cluster_j_time_g), main=paste0("Kernel ",1,", t=1") )


# ORA step=5 => n_time=320
phi_unlist <- unlist(phi)
num_plot <- 20
n2mfrow(num_plot)

{
  X11(width=1300, height=900)
  par(mfrow=n2mfrow(num_plot))
  par(oma=c(0,0,2,0))
  cluster <- 1
  for(time in 1:num_plot)
  {
    phi_cluster_1_time_time <- phi_unlist[seq(time + (cluster-1)*M , length(phi_unlist), by=M*n_time)]
    traceplot(as.mcmc(phi_cluster_1_time_time), 
              main=paste0("time = ",time))#, ylim=c(0,M))
  }
  title(paste0("phi_t^2"), outer = TRUE)
}


for(i in 1:(320/num_plot))
{
  X11(width=1300, height=900)
  par(mfrow=n2mfrow(num_plot))
  par(oma=c(0,0,2,0))
  cluster <- i
  for(time in 1:num_plot)
  {
    phi_cluster_1_time_time <- phi_unlist[seq(time + (cluster-1)*M , length(phi_unlist), by=M*n_time)]
    traceplot(as.mcmc(phi_cluster_1_time_time), 
              main=paste0( "time = ",time+20*(i-1) )
              )#, ylim=c(0,M))
  }
  title(paste0("phi_t^2"), outer = TRUE)
}

graphics.off()



## UNLIST ##--------------------------------------------------------
hey <- matrix(c(1,2,3,4,5,6), nrow=2, ncol=3)
yo <- hey-7
hey
yo

lis <- list(hey,yo)
lis

ul <- unlist(lis)
ul

length(ul) 

ul[seq(1,length(ul) ,by=3*2)] # by=30*150

# yeaaah





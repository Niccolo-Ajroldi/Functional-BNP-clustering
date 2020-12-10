
setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

rm(list=ls())
cat("\014")

library(invgamma)
library(fda)
library(MASS)

# M: number of clusters
# n_time: number of time points 
# n: number of observations
# L: number of basis
# n_iter: number of iterations

M <- 30
n <- 26
n_time <- 1600
L <- 60
n_iter <- 1

# basis<- functional data object FDA package


K <- matrix(0,n_iter,n_obs) 
#(K)_ij: cluster di assegnazione all'iterazione i dell'osservazione j


# Stick breaking
V <- numeric(M)
p <- rep(1/M, M) 
#vettore dei pesi inizializzato con un'uniforme

# total mass parameter
mass <- 2

# prior parameters of mu
m0 <- rep(0,L)
xi <- 1 # variance parameter da elicitare
Lamda0 <- diag(xi,L)
Lamda0_inv <- diag(1/xi,L) 

# prior parameters of sigma^2
a <- 1
b <- 1
  
# prior parameters of phi_t
c <- 1
d <- 1
  

#inializzazione parametri
mu <- matrix(M,n_time)
#(MU)_jt: valore di mu_j(t)
#i.e. media del cluster j valutata all'istante t (t=1:1600)
beta <- matrix(0,M,L)
#(beta)_ij: i-th kernel, j-th beta
sigma2 <- numeric(M)
#(sigma2)_i: sigma^2 del cluster i-esimo
phi <- matrix(0,M,n_time)
#(PHI_t)_ij: phi_t del cluster i-esimo al time point j-esimo (j=1:1600)

#dati
# X_obs <- matrix n_obs x n_time
#(X_obs)_ij: valore di X_i(j)
#i.e. osservazione i all'istante j (j=1:1600)

##  
L <- 60 # number of basis functions
m <- 4  # basis spline order (1+degree)
basis <- create.bspline.basis(rangeval=c(0,n_times), nbasis=L, norder=m)

time.grid <- 1:1600 # TODO: da riscalare
# b: matrix: L x n_times of basis fuctions evaliated on the time time grid
b <- t(eval.basis(time.grid, basis))
# b2: n_times vector, b2[t] is the dot product between b[,t] and itself
b2 <- sapply(1:n_times, function(t) crossprod(b[,t], b[,t]))

# smooth data
load("smooth_60b_nopenalization.RData")
X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)
# save coefficients
beta <- t(X_smoothed_f$fd$coefs)


n_iter <- 1
pb <- txtProgressBar(0, n_iter, style = 3)


for(iter in 1:n_iter)
{
  setTxtProgressBar(pb, iter)
  
  #### STEP 1 ####---------------------------------------------------------
  
  # (pi)j: probability that observation i belongs to cluster j
  p_i <- numeric(M)
  
  # iterate over the observations
  for(i in 1:n)
  {
    # iterate over kernels
    for(j in 1:M)
    {
      p_i[j] <- p[j] - 0.5*sum( log(2*pi*sigma2[j]*phi[j,]) + (X[i,]-mu[j,])^2/(2*sigma2[j]*phi[j,]) ) # Nan
    }
    p_i <- exp(p_i)/sum(exp(p_i))
    
    # cluster assignment at current iteration
    K[iter,i] <- sample(1:M, size=1, prob=p_i)
  }
  
  #### STEP 2 ####---------------------------------------------------------
  
  V[1] <- rbeta(1, 1 + sum(K[iter,] == 1), mass + sum(K[iter,] > 1))
  p[1] <- V[1]
  
  for(l in 2:(M - 1))
  {
    V[l] <- rbeta(1, 1 + sum(K[iter,] == l), mass + sum(K[iter,] > j))
    p[l] <- V[l] * prod(1 - V[1:(l - 1)])
  }
  
  V[M] <- 1
  p[M] <- V[M] * prod(1 - V[1:(M - 1)])
  
  
  #### STEP 3 ####---------------------------------------------------------
  
  # iterate over kernels
  for(j in 1:M)
  {
    # observations in cluster j
    indexes_j <- which(K[iter,]==j)
    # number of observations in cluster j
    r <- length(indexes_j)
    
    # if the kernel is not empty
    if(r>0)
    {
      ## SIGMA2
      a_r <- a + r*n_times*0.5
      b_r <- b
      for(g in indexes_j)
      {
        b_r <- b_r + 0.5*sum( (X[g,] - mu[j,])^2/phi[j,] )
      }
      sigma2[j] <- rinvgamma(n = 1, shape=a_r, scale=b_r)
      
      ## PHI
      c_r <- c + r*n_times*0.5
      for(t in 1:n_times)
      {
        d_r <- d + sum( (X[indexes_j,t] - mu[j,t])^2/(2*sigma2[j]) )
        phi[j,t] <- rinvgamma(n=1, shape=c_r, scale=d_r)
      }
      
      ## MU
      magic <- 1/sigma2 * sum( b2/phi[j,] ) * diag(L)
      
      Lambda_r <- solve(Lambda0_inv + r*magic) 
      mu_r <- solve(Lambda_r) %*% ( Lambda0_inv %*% m0 + magic %*% colSums(beta[indexes_j,]) ) # funziona
    
      # sample coefficients of basis projection
      mu_coef <- mvrnorm(n=1, mu=mu_r, Sigma=Lambda_r)
      
      # evaluation of mu on the time grid
      mu[j,] <- mu_coef %*% b
      
    }
    
    else
    {
      ## SIGMA2
      sigma2[j] <- rinvgamma(n = 1, shape=a, scale=b)
      
      ## PHI
      phi[j,] <- rinvgamma(n=n_times, shape=c, scale=d)
      
      ## MU
      # sample coefficients of basis projection
      mu_coef <- mvrnorm(n=1, mu=m0, Sigma=Lamda0)
      # evaluation of mu on the time grid
      mu[j,] <- mu_coef %*% b
      
    }
    
  }
  

}





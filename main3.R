
#setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
#setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd('C:/Users/Utente/Desktop/Bayesian statistics/Functional-BNP-clustering')
#setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

rm(list=ls())
cat("\014")

library(invgamma)
library(fda)
library(MASS)
library("tictoc")

# M: number of clusters
# n_time: number of time points 
# n: number of observations
# L: number of basis
# n_iter: number of iterations

M <- 150
n <- 22
n_time <- 1600
L <- 40
n_iter <- 10000

#(K)_ij: cluster di assegnazione all'iterazione i dell'osservazione j
K <- matrix(0,nrow=n_iter,ncol=n) 

# Stick breaking
V <- numeric(M)
p <- rep(1/M, M) #vettore dei pesi inizializzato con un'uniforme

# total mass parameter
mass <- 10

# prior parameters of mu
m0 <- rep(0,L)
xi <- 1000 # variance parameter da elicitare
Lambda0 <- diag(xi,L)
Lambda0_inv <- diag(1/xi,L) 

# prior parameters of sigma^2
a <- 3
b <- 100

# prior parameters of phi_t
c <- 3
d <- 100

#inializzazione parametri
mu <- matrix(0,M,n_time)
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
m <- 4  # basis spline order (1+degree)
basis <- create.bspline.basis(rangeval=c(0,n_time), nbasis=L, norder=m)

time.grid <- 1:1600 # TODO: da riscalare
# basis.t: matrix: L x n_time of basis fuctions evaliated on the time time grid
basis.t <- t(eval.basis(time.grid, basis))
# basis2.t: n_time vector, b2[t] is the dot product between basis.t[,t] and itself
basis2.t <- sapply(1:n_time, function(t) crossprod(basis.t[,t], basis.t[,t]))

# smooth data
load("smooth_60b_nopenalization.RData")
 X <- X[-c(12,13,19),]
 X<-X[-21,]
X <- X/max(X) # TODO: capire se si può migliorare
X_smoothed_f <- smooth.basis(argvals=time.grid, y=t(X), fdParobj=basis)
# save coefficients
beta <- t(X_smoothed_f$fd$coefs)

#### Initialization ####

K[1,] <- sample(1:n)

## SIGMA2
sigma2 <- rinvgamma(n = M, shape=a, scale=b)

## PHI
phi <- matrix(rinvgamma(n=M*n_time, shape=c, scale=d), byrow=TRUE, nrow=M, ncol=n_time)

## MU
# sample coefficients of basis projection
mu_coef <- matrix(mvrnorm(n=M, mu=m0, Sigma=Lambda0), byrow=TRUE, nrow=M, ncol=L)
# evaluation of mu on the time grid
mu <- matrix(mu_coef %*% basis.t, byrow=TRUE, nrow=M, ncol=n_time)


# logarithmic scale translation
#C <- -1e7

pb <- txtProgressBar(0, n_iter, style = 3)
#pb <- txtProgressBar(0, n_iter, style = "ETA")

tic()

for(iter in 1:n_iter)
{
  
  #### STEP 1: acceleration step ####--------------------------------------------
  
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
      sigma2_j <- sigma2[j]
      phi_j <- phi[j,]
      mu_j <- mu[j,]
      
      ## SIGMA2
      a_r <- a + r*n_time*0.5
      b_r <- b
      for(g in indexes_j)
      {
        b_r <- b_r + 0.5*sum( (X[g,] - mu[j,])^2/phi[j,] )
      }
      sigma2[j] <- rinvgamma(n = 1, shape=a_r, scale=b_r)
      
      ## PHI
      c_r <- c + r*n_time*0.5
      for(t in 1:n_time)
      {
        d_r <- d + sum( (X[indexes_j,t] - mu[j,t])^2/(2*sigma2_j) )
        phi[j,t] <- rinvgamma(n=1, shape=c_r, scale=d_r)
      }
      
      ## MU
      magic <- 1/sigma2_j * sum( basis2.t/phi_j )
      
      Lambda_r <- solve(Lambda0_inv + r*magic*diag(L))
      if(r>1)
      {
        mu_r <- Lambda_r %*% ( Lambda0_inv %*% m0 + magic * colSums(beta[indexes_j,]) ) # funziona
      } else {
        mu_r <- Lambda_r %*% ( Lambda0_inv %*% m0 + magic * beta[indexes_j,] )
      }
      
      # sample coefficients of basis projection
      mu_coef <- mvrnorm(n=1, mu=mu_r, Sigma=Lambda_r)
      
      # evaluation of mu on the time grid
      mu[j,] <- mu_coef %*% basis.t
      
    } else {
      ## SIGMA2
      sigma2[j] <- rinvgamma(n = 1, shape=a, scale=b)
      
      ## PHI
      phi[j,] <- rinvgamma(n=n_time, shape=c, scale=d)
      
      ## MU
      # sample coefficients of basis projection
      mu_coef <- mvrnorm(n=1, mu=m0, Sigma=Lambda0)
      # evaluation of mu on the time grid
      mu[j,] <- mu_coef %*% basis.t
      
    }
    
  }
  
  #### STEP 2: K ####---------------------------------------------------------
  
  # (pi)j: probability that observation i belongs to cluster j
  p_i <- numeric(M)
  
  # iterate over the observations
  for(i in 1:n)
  {
    # iterate over kernels
    for(j in 1:M)
    {
      p_i[j] <- log(p[j]) + sum( (-0.5)*log(2*pi*sigma2[j]*phi[j,]) + (X[i,]-mu[j,])^2/(2*sigma2[j]*phi[j,]) ) # Nan
    }
    p_i <- p_i - max(p_i) # TODO: controllare
    p_i <- exp(p_i)/sum(exp(p_i))
    # cluster assignment at current iteration
    K[iter,i] <- sample(1:M, size=1, prob=p_i)
  }
  
  #### STEP 3: weights ####---------------------------------------------------------
  
  V[1] <- rbeta(1, 1 + sum(K[iter,] == 1), mass + sum(K[iter,] > 1))
  p[1] <- V[1]
  
  for(l in 2:(M - 1))
  {
    V[l] <- rbeta(1, 1 + sum(K[iter,] == l), mass + sum(K[iter,] > j))
    p[l] <- V[l] * prod(1 - V[1:(l - 1)])
  }
  
  V[M] <- 1
  p[M] <- V[M] * prod(1 - V[1:(M - 1)])
  
  ##---------------------------------------------------------------------------------
  
  
  setTxtProgressBar(pb, iter)
  
}

toc()

View(K)

save(K, file=paste0("K_", n_iter, "_iter_no21_b",b,".RData"))

ind <- 4993
row <- K[ind,]
unique(row)

x11()
matplot(1:1600,t(X), type='l', xlab='x', ylab='orig.func', col=row)

s <- apply(K, 1, function(x) length(unique(x))>1)
su <- unlist(s)
rows <- which(s==TRUE)

r <- sample(rows,10)
for(i in r)
{
  x11()
  matplot(1:1600,t(X), type='l', xlab='x', ylab='orig.func', col=K[i,])
  legend("topright", legend = levels(factor(K[i,])) , col=unique(K[i,]))#, col=unique(K[i,]),cex=0.8)
}

graphics.off()





# utility

library(invgamma)

m <- 100
v <- 10
a <- m^2/v + 2
b <- m^3/v + m

a <- c(0.01,0.05,0.1,0.5,1,5,10)
b <- 0.1
x11()
par(mfrow=c(3,3))
for(i in 1:length(a))
{
  y <- rgamma(100000, a[i], b)
  plot(density(y), col="lightblue", lwd=3, main=paste0("Shape parameter: ",a[i]))
}
abline(v=m, lwd=1, col="red")


#### Nuove osservazioni --> 3 cluster meno separati
mu <- mu[j,]
X <- X[i,]

#phi[j,] <- rinvgamma(n=n_time, shape=c, rate=d)

fun <- function(PHI)
{
  #val <- sum((-0.5)*log(2*pi*phi[j,]) - ((X[i,]-mu[j,])^2)/(2*phi[j,]) )
  val <- sum((-0.5)*log(2*pi*PHI)*rep(1:30) - ((X-mu)^2)/(2*PHI) )
  return(val)
}


PHI <- seq(0.0001,10,by=0.0001)
y <- numeric(length(PHI))
for(i in 1:length(PHI))
{
  y[i] <- fun(PHI[i])
}

x11()
plot(PHI,y, type="l", lwd=3, col="lightblue")

# utility

m <- 100
v <- 10
a <- m^2/v + 2
b <- m^3/v + m
y <- rinvgamma(100000, a, b)
x11()
plot(density(y), col="lightblue", lwd=3)
abline(v=m, lwd=1, col="red")


#### Nuove osservazioni --> 3 cluster meno separati
mu <- mu[j,]
X <- X[i,]

#phi[j,] <- rinvgamma(n=n_time, shape=c, rate=d)

fun <- function(PHI)
{
  #val <- sum((-0.5)*log(2*pi*phi[j,]) - ((X[i,]-mu[j,])^2)/(2*phi[j,]) )
  val <- sum((-0.5)*log(2*pi*PHI)*rep(1:100) - ((X-mu)^2)/(2*PHI) )
  return(val)
}


PHI <- seq(0.001,10,by=0.001)
y <- numeric(length(PHI))
for(i in 1:length(PHI))
{
  y[i] <- fun(PHI[i])
}

x11()
plot(PHI,y, type="l", lwd=3, col="lightblue")

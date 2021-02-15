
# number of observations
n <- dim(X)[1]
  
# total mass parameter
mass <- 0.65

# expected number of clusters
sum(mass/(mass+0:(n-1)))

# approximation 1
mass*log(1+n/mass)

# approximation 2
mass*log(n)


# Approximate distance between the marginal density of X1,...,Xn with the truncated stick breaking
# and the marginal density of X1,...,Xn with the infinite stick breaking
# [Ishwaran & James theorem 2]
M <- 50
mass <- 2
n <- 26
4*n*exp(-(M-1)/mass) 
# NOPEEEEEEEEEEEEEE (infinite dimensional problem, I&W è pe runidimansionale, no?)


# secondo me non ha più molto senso mettere massa tale per cui vengano solo 2 custer
# perchè abbiam notato che anche i GOSE son molto diversi uno dall'altro
# ocio però che se aumentiamo la massa, anche il numero di cluster deve aumentare, per avere una buona approsimazione
# quindi allora ha senso tenere alfa non troppo alto (in modo da avere 5/6 cluster attesi)
# in modo da poter tenere basso anche M e avere un algoritmo più veloce!

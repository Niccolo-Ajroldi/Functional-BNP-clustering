
# number of observations
n <- dim(X)[1]
  
# total mass parameter
mass <- 0.31

# expected number of clusters
sum( mass/(mass+(0:(n-1))) )

# approximation 1
mass*log(1+n/mass)

# approximation 2
mass*log(n)

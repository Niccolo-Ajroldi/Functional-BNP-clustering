
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

#load("Results/Nico_31_12.RData")

source("fun_traceplots.R")

# save old directory to get back there at the end
dir.current <- getwd()

fun_traceplots (out,
                smoothing_list,
                run_parameters,
                blocks=2,
                time=0.15,
                mu_gif=TRUE,
                iter_step=20)

setwd(dir.current)


clust <- 3

K <- out$K
View(K)

clust <- c(1,3,4,5,6,9,10,17,21,22)

x11()
par(mfrow=c(2,1))
matplot(t(X[clust,]), type="l")
matplot(t(X[-clust,]), type="l")


mu_coef_out <- out$mu_coef_out
is(mu_coef_out)
dim(mu_coef_out[[1]])
mu_coef_out[[1]][1,]
mu_coef_out[[2]][1,]
mu_coef_out[[3]][1,]
mu_coef_out[[4]][1,]
mu_coef_out[[5]][1,]
mu_coef_out[[6]][1,]
mu_coef_out[[7]][1,]
mu_coef_out[[8]][1,]
mu_coef_out[[9]][1,]
mu_coef_out[[10]][1,]


for(i in 1:50)
{
  x11()
  plot(mu_coef_out[[i]][1,])
}

graphics.off()

time.grid <- smoothing_list$time.grid
basis <- smoothing_list$basis
basis.t <- t(eval.basis(time.grid, basis))

for(i in 1:50)
{
  x11()
  matplot(t(mu_coef_out[[i]][1,] %*% basis.t), type="l")
}
graphics.off()

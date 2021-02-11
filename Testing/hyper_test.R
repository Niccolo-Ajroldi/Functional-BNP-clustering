library(invgamma)

var_phi.grid <- c(0.0001, 0.001, 0.01, 1,
                 10, 100, 1000, 10000)
for(var_phi in 1){
  hyper_list <- hyperparameters(var_phi = var_phi, 
                                X = smoothing_list$X,
                                beta = smoothing_list$beta,
                                scale = 1)
  
  plot(seq(0,50,by=0.01), dinvgamma(seq(0,50,by=0.01),
                                    shape=hyper_list$c,
                                    rate=hyper_list$d),
       type='l',
       main=paste0('Varphi=',var_phi))
  abline(v=10, col='red')
}



sigma.var.grid <- seq(0.1,2, length=50)
#phi.var.grid <- seq(0.1,10, length=10)

for(jj in 1:length(sigma.var.grid))
{
  print(jj)
  sigma.var <- sigma.var.grid[jj]
  phi.var   <- sigma.var.grid[jj]
  out <- FBNP(n_iter = 20,
              burnin = 0,
              thin = 1,
              M = 1000,
              mass = 2,
              smoothing = smoothing_list,
              hyperparam = hyper_list)
  run_parameters <- list('algorithm_parameters' = out$algorithm_parameters,
                         'prior_parameters'     = hyper_list,
                         'smoothing_parameters' = smoothing_list$smoothing_parameters)
  traceplot_K(out, 
              smoothing_list, 
              run_parameters,
              plot.title = paste0("Var(sigma)=",sigma.var)
              )  
}

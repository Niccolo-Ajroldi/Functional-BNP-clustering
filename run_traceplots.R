
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")

source("fun_traceplots.R")

# save old directory to get back there at the end
dir.current <- getwd()

fun_traceplots (out,
                smoothing_list,
                run_parameters,
                blocks=10,
                time=0.15,
                mu_gif=TRUE,
                iter_step=10)

setwd(dir.current)

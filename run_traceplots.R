
# setwd("C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/Functional-BNP-clustering')
setwd("D:/Poli/Corsi/BAYESIAN/Proj/Functional-BNP-clustering")
# setwd('C:/Users/edoar/Desktop/Bayesian statistics/Project/code/No github code')

#load("Results/Nico_31_12.RData")

source("fun_traceplots.R")

# save old directory to get back there at the end
dir.current <- getwd()

# call
fun_traceplots (out,
                smoothing_list,
                run_parameters,
                blocks=1,    # each "block" contains 10 kernels, change this in order to see a different gif for mu(t)
                time=0.15,   # time between images in the gif
                mu_gif=FALSE,
                falzo=FALSE  # change this if you are using falzo (different scale of mu(t))
                )

# go to previous directory
setwd(dir.current)



graphics.off()

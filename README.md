# Functional-BNP-Clustering
###  Bayesian nonparametric functional clustering

This is the repository for the project of the course *Bayesian Statistics* held by Professor A. Guglielmi at Politecnico di Milano during the academic year 2020/2021.
The tutor of the project was R. Corradin, PhD.

We were given 26 multivariate functional observation with 4 components, each of which represents a somatosensory evoked potential, i.e. the electrical modifications occurring in the central nervous system following a stimulus, detected within 2 seconds after the stimulus over 1600 instants.

The aim of the project is **to implement a functional clustering algorithm in a Bayesian nonparametric framework** in a univariate setting. In particular, we assume an infinite mixture model for our functional obervations with a Dirichlet Process as mixture distribution. Our aim is to sample from Dirichlet process in order to get a sample of the latent partition induced. A Gibbs sampler is implemented for our purposes.

The algorithm and the model will be first tested on simulated data and then on clinical ones.

A short explanation of the model is made in this file, for further details and for a deep mathematical explanation please refer to the [report](link).

The core of the Gibbs Sampler is implemented in the function `FBNP_hyper.R`.

### Model in a nutshell

Our observed functions are assume to be realizations of random functions that follow a mixture of Gaussian Processes, with a Dirichlet Process as mixture distribution. We are interested in the latent random parameters defyining the Gaussian Process of each observation and we will exploit ties to cluster them. We perform a dimensionality reduction assuming stationarity and independence between different time points. On the other hand, the mean operator is expanded over a basis that we assume exactly spans the space of our data.  

### Sampling algorithm 

Our sampling algorithm is a Blocked Gibbs sampler presented in 
> Ishwaran, H., and James, L.F. (2001). “Gibbs sampling methods for stick-breaking priors”. Journal of the American Statistical Association, 96.453, pp. 161–173.

It is noteworthy that since we consider an hyperprior structure, we will refer to an slightly modified version of the Blocked Gibbs sampler algorithm, taking into account the full conditional for the hyperparameters.


## Repository structure
The repository is presents the following elements:
* `main.R`: R-script with all the steps for running the Blocked Gibbs sampler
* `Tools`: folder containing all the script needed in `main.R`, among whom `FBNP.R` containing a function that implements Gibbs sampler **without hyperpriors** and
     `FBNP_hyper.R` containing a function that implements Gibbs Sampler **with hyperpriors**
* `Posterior inference`: folder contanining the script `main_posterior inference.R` for performing posterior inference for the latent partition and convergence diagnostic and all the script neeeded in the main one
* `Gaussian process`: folder containing the scripts `main_GP_indep.R` and `main_GP_exp.R` that perform the Blocked Gibbs sampler on the simulated data with diagonal and exponential covariance operator and all the needed script
* `Alltime`: folder containing `main_alltime.R` that performs the Blocked Gibbs sampler changing the prior structure on the covariance operator, setting a different prior on different time instants
* `Ajroldi-Bortolotti-Marchionni_Functional-BNP-Clustering.pdf`: report of the project

### Installation

The file `install.R` provides automatic installation of the required packages.

### Simulated data
We simulate 10 trajectories of 3 different Gaussian Processes with different sinusoidal mean operator and a covariance operator coherent with our model (for further details please refer to [report](link)).
After having conveniently tuned the parameters, our algorithm is able to properly cluster the simulated data with the hyperprior structure.
![alt text](https://github.com/Niccolo-Ajroldi/Functional-BNP-clustering/blob/main/pics/GP_ind.png)




### Clinical data
Since we work in a univariate setting, we select one of the four components of our functional observations. In particular, we run our algorithm on the right lobe short latency signal and then for the other components
![alt text](https://github.com/Niccolo-Ajroldi/Functional-BNP-clustering/blob/main/pics/Data_cutted.png)
For posterior inference on the latent partition please refer to the  [report](link).

### References
Computation were performed using 
>  R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

For all the R-packages used in the codes please refer to section Reference in the [report](link).



## Authors

Niccolò Ajroldi - Politecnico di Milano  
Teresa Bortolotti - Politecnico di Milano 
Edoardo Marchionni - Politecnico di Milano
(Master of Science in Mathematical engineering students at Politecnico di Milano)


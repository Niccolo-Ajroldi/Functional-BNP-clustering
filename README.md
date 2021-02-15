# Functional-BNP-Clustering
###  Bayesian nonparametric functional clustering

This is the repository for the project of the course *Bayesian Statistics* held by Professor A. Guglielmi at Politecnico di Milano during the academic year 2020/2021.
The tutor of the project was R. Corradin, PhD.

We were given 26 multivariate functional observation with 4 components, each of which represents a somatosensory evoked potential, i.e. the electrical modifications occurringin the central nervous system following a stimulus, detected within 2 seconds after the stimulus over 1600 instants.

The aim of the project is **to implement a functional clustering algorithm in a Bayesian nonparametric framework** in a univariate setting. In particular, we assume an infinite mixture model for our functional obervations with a Dirichlet Process as mixture distribution. Our aim is to sample from Dirichlet process in order to get a sample of the latent partition induced. A Gibbs sampler is implemented for our purposes.

The algorithm and two different version of model, with and without hyperpriors for the parameters of the base measure of the Dirichlet Process,  will be first tested on simulated data and then on clinical ones.

A short explanation of the model is made in this file, for further details and for a deep mathematical explanation please refer to the [report](link).

The core of the Gibbs Sampler is implemented in the function `FBNP.R`.

### Model in a nutshell

Our observed functions are assume to be realizations of random functions that follow a mixture of Gaussian Processes, with a Dirichlet Process as mixture distribution. We are interested in the latent random parameters defyining the Gaussian Process of each observation and we will exploit ties to cluster them. We perform a dimensionality reduction assuming stationarity and independence between different time points. On the other hand, the mean operator is expanded over a basis that we assume exactly spans the space of our data.  

### Sampling algorithm 

Our sampling algorithm is a Blocked Gibbs sampler presented in 
> Ishwaran, H., and James, L.F. (2001). “Gibbs sampling methods for stick-breaking priors”. Journal of the American Statistical Association, 96.453, pp. 161–173.

It is noteworthy that in the case with hyperpriors, we will refer to an slightly modified version of the Blocked Gibbs sampler algorithm, taking into account the full conditional for the hyperparameters.


## Repository structure
The repository is presents the following elements:
* `main.R`: R-script with all the steps for clustering ???
* `FBNP`: folder containing script 
    * `FBNP.R` containing a function that implements Gibbs Sampler **without hyperpriors**
    * `FBNP_hyper.R` containing a function that implements Gibbs Sampler **with hyperpriors**
* `Smoothing`
* `Hyperparameters`
* `Ajroldi-Bortolotti-Marchionni_Functional-BNP-Clustering.pdf`: report of the project

### Installation

The file `install.R` provides automatic installation of the required packages, which are reported in Section??.

### Simulated data

![alt text](https://github.com/Niccolo-Ajroldi/Functional-BNP-clustering/blob/main/pics/Simulated_GP.png)


### Clinical data
![alt text](https://github.com/Niccolo-Ajroldi/Functional-BNP-clustering/blob/main/pics/Data_cutted.png)

### References
Computation were performed using 
>  R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

For all the R-packages used in the codes please refer to section Reference in the [report](link).

This markdown file is written relying on [DILLINGER](https://dillinger.io/).

## Authors

Niccolò Ajroldi - Politecnico di Milano  
Teresa Bortolotti - Politecnico di Milano 
Edoardo Marchionni - Politecnico di Milano
(Master of Science in Mathematical engineering students at Politecnico di Milano)


# Functional-BNP-clustering

In this repository we implement a Dirichlet Process Mixture model to perform Bayesian Nonparametric clustering  of Functional Data. A detailed mathematical explanation of the model is addressed in the dedicated [report](). \
The core of the Gibbs Sampler is implemented in the function `FBNP.R`.

### Model

Observed functions x1(t),...,xn(t) are assume to be realizations of random functions X1(t),...,Xn(t). The probabiliy distribution of such random functions is a mixture of Gaussian Processes, with Dirichlet Process as mixing measure. \
We are interested in the latent random variables defyining each GP. We will exploit ties between them to cluster observations. \
We report here the model used for the algorithm, for a more detailed explanation refer to the attached report.

### Repository structure

Folders and codes are structured in the following way:
- smoothing: codes for smoothing of functional data
- GP_simulation: codes for testing the model on data simulated from Gaussian Processes
- ...
- `FBNP.R` implements the Gibbs Sampler **without hyperpriors**
- `FBNP_hyper.R` implements the Gibbs Sampler **with hyperprior**
- Ajroldi-Bortolotti-Marchionni-Functional-BNP-Clustering.pdf: contains a report of the project

### Installation

The file `install.R` provides automatic installation of the required packages, which are reported in Section??.

### Simulated data

![alt text](https://github.com/Niccolo-Ajroldi/Functional-BNP-clustering/blob/main/pics/Simulated_GP.png)


### Clinical data

### References

Codes are written in [R](https://www.r-project.org/): A language and environment for statistical computing. \
The following packages have been used in this project:
- [LaplaceDemon](https://web.archive.org/web/20150206004624/http://www.bayesian-inference.com/software) Statisticat, LLC. (2020). LaplacesDemon: Complete Environment for Bayesian Inference. R package version 16.1.4.
- [coda](https://journal.r-project.org/archive/) Plummer M, Best N, Cowles K, Vines K (2006). “CODA: Convergence Diagnosis and Output Analysis for MCMC.” R News, 6(1), 7–11.
- [roahd](https://CRAN.R-project.org/package=roahd) Ieva F, Paganoni AM, Romo J, Tarabelloni N (2019). “roahd Package: Robust Analysis of High Dimensional Data.” The R Journal, 11(2), 291–307.
- [MASS](http://www.stats.ox.ac.uk/pub/MASS4/) Venables WN, Ripley BD (2002). Modern Applied Statistics with S, Fourth edition. Springer, New York. ISBN 0-387-95457-0.
- [pbmcapply](https://CRAN.R-project.org/package=pbmcapply )
- [fda](https://cran.r-project.org/web/packages/fda/index.html)
- [fdakma](https://cran.r-project.org/web/packages/fdakma/index.html)
- [invgamma](https://CRAN.R-project.org/package=invgamma)

This markdown file is written relying on [DILLINGER](https://dillinger.io/).

### Authors

Niccolò Ajroldi - Politecnico di Milano \
Teresa Bortolotti - Politecnico di Milano \
Edoardo Marchionni - Politecnico di Milano


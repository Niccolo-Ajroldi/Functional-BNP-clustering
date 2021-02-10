# Functional-BNP-clustering

In this repository we implement a Dirichlet Process Mixture model to perform Bayesian Nonparametric clustering  of Functional Data. \
A detailed mathematical explanation of the model is addressed in the [report](). \
The core of the Gibbs Sampler is implemented in the function `FBNP.R`.

### Model

... Spiegazione modello ...

Qui mettiamo una bella foto del modello presa da Latex oppure riportiamo il codice latex in markdown.

### Repository structure

Folders and codes are structured in the following way:
- smoothing: codes for smoothing of functional data
- GP_simulation: codes for testing the model on data simulated from Gaussian Processes
- ...
- Ajroldi-Bortolotti-Marchionni-Functional-BNP-Clustering.pdf: contains a report of the project

### Installation

The file `install.R` provides automatic installation of the required packages, which are reported in Section??.

### Simulated data

### Clinical data

### References

Codes are written in [R](https://www.r-project.org/): A language and environment for statistical computing. \
The following packages have been used in this project:
- [LaplaceDemon](https://CRAN.R-project.org/package=LaplacesDemon) Statisticat, LLC. (2020). LaplacesDemon: Complete Environment for Bayesian Inference. R package version 16.1.4, https://web.archive.org/web/20150206004624/http://www.bayesian-inference.com/software.
- [coda](https://CRAN.R-project.org/package=coda) Plummer M, Best N, Cowles K, Vines K (2006). “CODA: Convergence Diagnosis and Output Analysis for MCMC.” R News, 6(1), 7–11. https://journal.r-project.org/archive/.
- [fda](https://cran.r-project.org/web/packages/fda/index.html)
- [roahd](https://CRAN.R-project.org/package=roahd) Ieva F, Paganoni AM, Romo J, Tarabelloni N (2019). “roahd Package: Robust Analysis of High Dimensional Data.” The R Journal, 11(2), 291–307. https://doi.org/10.32614/RJ-2019-032
- [fdakma]()
- [invgamma](https://CRAN.R-project.org/package=invgamma)
- [MASS](https://CRAN.R-project.org/package=MASS) Venables WN, Ripley BD (2002). Modern Applied Statistics with S, Fourth edition. Springer, New York. ISBN 0-387-95457-0 http://www.stats.ox.ac.uk/pub/MASS4/.
- [pbmcapply](https://CRAN.R-project.org/package=pbmcapply )

This markdown file is written relying on [DILLINGER](https://dillinger.io/).

### Authors

Niccolò Ajroldi - Politecnico di Milano \
Teresa Bortolotti - Politecnico di Milano \
Edoardo Marchionni - Politecnico di Milano


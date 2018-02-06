QUBIC
=====

[![Build Status](https://travis-ci.org/zy26/QUBIC.svg?branch=devel)](https://travis-ci.org/zy26/QUBIC)

Authors: [Yu Zhang](mailto:zy26@jlu.edu.cn) and [Qin Ma](mailto:qin.ma@sdstate.edu)

Desription
------------
QUBIC is recognized as one of the best biclustering methods in terms of its efficiency and effetiveness in biological data interpretion. This package provides an R implementation of the QUBIC algorithm, with significantly improved efficiency and comprehensive functions. 

Installation
------------

Please follow the instructions mentioned in the URL: http://bioconductor.org/packages/QUBIC

To install the development version of QUBIC, you will need to install at least the following packages from CRAN
```{r}
install.packages("biclust")
install.packages("Rcpp")
install.packages("RcppArmadillo")
source("http://bioconductor.org/biocLite.R") # install BiocInstaller
```
For Windows users, Rtools(https://cran.r-project.org/bin/windows/Rtools/) should also be installed.

Then,
```{r}
install.packages("devtools")
devtools::install_github("zy26/QUBIC")
```

To draw heatmap and visualize network, the following packages should also be installed:
```{r}
install.packages("qgraph")
install.packages("RcolorBrewer")
```

Vignette
------------
You can find the tutorial for QUBIC at https://bioconductor.org/packages/release/bioc/vignettes/QUBIC/inst/doc/qubic_vignette.pdf

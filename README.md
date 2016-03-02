QUBIC
=====

[![Build Status](https://travis-ci.org/zy26/QUBIC.svg?branch=master)](https://travis-ci.org/zy26/QUBIC)

Authors: [Yu Zhang](mailto:zy26@jlu.edu.cn) and [Qin Ma](mailto:qin.ma@sdstate.edu)

Type: R package


Desription
------------
QUBIC is recognized as one of the best biclustering methods in terms of its efficiency and effetiveness in biological data interpretion. This package provides an R implementation of the QUBIC algorithm, with significantly improved efficiency and comprehensive functions. 

Installation
------------
First, you will need to install at least the following packages from CRAN:
```{r}
install.packages("biclust")
install.packages("Rcpp")
install.packgaes("qgraph")
install.packages("RColorBrewer")
```
For Windows users, Rtools(https://cran.r-project.org/bin/windows/Rtools/)should also be installed.

To install the development version of QUBIC:
```{r}
install.packages("devtools")

devtools::install_github("zy26/QUBIC")
```

Vignette
------------
You can find the vignette for QUBIC at https://github.com/zy26/QUBIC/blob/master/vignettes/qubic_vignette.Rmd

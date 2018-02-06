QUBIC
=====

Authors: [Yu Zhang](mailto:zy26@jlu.edu.cn) and [Qin Ma](mailto:qin.ma@sdstate.edu)

Description
------------

This package provides our R/Bioconductor implementation of the original QUalitative BIClustering (QUBIC) algorithm (Li et al., Nucleic Acids Research, [10.1093/nar/gkp491](https://doi.org/10.1093/nar/gkp491)), along with utilities for discretization, heatmap visualization, and bicluster-based co-expression network analysis.

In comparative benchmark studies under specific biclustering settings, QUBIC showed strong performance, including high enriched bicluster ratios on real datasets and robust behavior under increasing noise (Eren et al., Briefings in Bioinformatics, [10.1093/bib/bbs032](https://doi.org/10.1093/bib/bbs032)).

**Citing us:** Yu Zhang, Juan Xie, Jinyu Yang, Anne Fennell, Chi Zhang, Qin Ma. QUBIC: a Bioconductor package for qualitative biclustering analysis of gene co-expression data. Bioinformatics. 2017;33(3):450-452. <https://doi.org/10.1093/bioinformatics/btw635>

Installation
------------

Stable Bioconductor release:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager")
}
BiocManager::install("QUBIC")
```

Development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("zy26/QUBIC")
```

Optional visualization packages:

```r
install.packages("qgraph")
install.packages("RColorBrewer")
```

Notes:

- On Windows, install Rtools when building from source: <https://cran.r-project.org/bin/windows/Rtools/>
- `biclust` is optional for compatibility workflows and can be installed from CRAN Archive if needed.

Vignette
------------

You can find the tutorial for QUBIC at:
<https://bioconductor.org/packages/release/bioc/vignettes/QUBIC/inst/doc/qubic_vignette.pdf>

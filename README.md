# The Promise of Dynamic Phylogenetic Models for Cultural Evolution

Slides, abstract, and R code for talk at the Max Planck Institute for 
Evolutionary Anthropology on 26th November 2025.

To run the R code in `mpi-eva-talk.R`, you will need to 
[install R](https://www.r-project.org/) and the following R packages:

```r
install.packages(c("ape", "BiocManager", "brms", "cmdstanr",
                   "devtools", "tidybayes", "tidyverse"))
devtools::install_github(c("ScottClaessens/coevolve", "rgriff23/btw"))
BiocManager::install("ggtree")
```

You will also need to download `BayesTraitsV3.exe` into this repository
(see [here](https://www.evolution.reading.ac.uk/BayesTraitsV3/BayesTraitsV3.html)).

## Authors

Scott Claessens, scott.claessens@gmail.com

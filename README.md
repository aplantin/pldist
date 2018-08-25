# pldist: Paired and Longitudinal Ecological Dissimilarities 
**Author:** Anna Plantinga

## Introduction

`pldist` allows distance-based analysis of paired and longitudinal microbiome data. In particular, the package supports both paired and longitudinal versions of unweighted UniFrac, generalized UniFrac, Bray-Curtis, Jaccard, Gower, and Kulczynski distances or dissimilarities. Functions implementing the transformations that underlie these distances are also provided so that transformed OTU data may be included in analyses beyond distance-based methods. The code can handle paired data, balanced longitudinal data, and unbalanced longitudinal data, although use for highly unbalanced designs is not recommended.

## Issues

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/aplantin/pldist/issues). 

## Installation 

You may install `pldist` from GitHub using the following code: 

```{r install} 
## install.packages("devtools") # only run this line if necessary
devtools::install_github(repo = "aplantin/pldist")
```

## Example

This example demonstrates usage of `pldist` in a simple setting with simulated data. For more examples and details, please see the [vignette](https://github.com/aplantin/pldist/blob/master/vignettes/introduction.Rmd). 

```{r example} 
library(pldist)
library(ape) 
data(sim.tree)
data(paired.otus)
data(paired.meta)
``` 


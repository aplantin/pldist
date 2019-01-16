# pldist: Paired and Longitudinal Ecological Dissimilarities 
**Author:** Anna Plantinga

## Introduction

`pldist` allows distance-based analysis of paired and longitudinal microbiome data. In particular, the package supports both paired and longitudinal versions of unweighted UniFrac, generalized UniFrac, Bray-Curtis, Jaccard, Gower, and Kulczynski distances or dissimilarities. Functions implementing the transformations that underlie these distances are also provided so that transformed OTU data may be included in analyses beyond distance-based methods. The code can handle paired data, balanced longitudinal data, and unbalanced longitudinal data, although use for highly unbalanced designs is not recommended.

## Issues

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/aplantin/pldist/issues). 

## Installation 

You may install `pldist` from GitHub using the following code: 

```{r install} 
## install.packages("devtools") # run this line if devtools is not installed
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
data(bal.long.otus) 
data(bal.long.meta) 

# Look at the OTU data: 
# row names are sample IDs, column names are OTU IDs 
paired.otus[1:4,1:4] 

# Look at the metadata: 
# columns are subjID, sampID, time 
# One row per sample 
paired.meta[1:4, ]

# Gower distance, paired quantitative transformation, with/without CLR 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, clr = FALSE, method = "gower")$D
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, clr = TRUE, method = "gower")$D

# Gower distance, paired qualitative/binary transformation 
pldist(paired.otus, paired.meta, paired = TRUE, binary = TRUE, method = "gower")$D

# Gower distance, longitudinal quantitative transformation 
pldist(bal.long.otus, bal.long.meta, paired = FALSE, binary = FALSE, method = "gower")$D

# Gower distance, longitudinal qualitative/binary transformation 
pldist(bal.long.otus, bal.long.meta, paired = FALSE, binary = TRUE, method = "gower")$D

# Other distances 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "bray")$D
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "kulczynski")$D
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "jaccard")$D

# UniFrac also requires a phylogenetic tree and gamma values 
# (Gamma controls weight placed on abundant lineages) 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, 
    method = "unifrac", tree = sim.tree, gam = c(0, 0.5, 1))$D 
    
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, clr = TRUE, 
    method = "unifrac", tree = sim.tree, gam = c(0, 0.5, 1))$D 
``` 


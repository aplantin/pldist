# pldist: Paired and Longitudinal Ecological Dissimilarities 
**Author:** Anna Plantinga

# Introduction

# Issues

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/aplantin/pldist/issues). 

# Installation 

You may install `pldist` from GitHub using the following code: 

```{r install} 
## install.packages("devtools") # only run this line if necessary
devtools::install_github(repo = "aplantin/pldist")
```

# Example

This example demonstrates usage of `pldist` in a simple setting with simulated data. For more examples and details, please see the [vignette](https://github.com/aplantin/pldist/blob/master/vignettes/introduction.Rmd). 

```{r example} 
library(pldist)
library(ape) 
data(sim.tree)
data(paired.otus)
data(paired.meta)
``` 


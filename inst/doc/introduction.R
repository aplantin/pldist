## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installation-instructions-------------------------------------------
# only run if you don't have devtools installed 
#install.packages("devtools"); library(devtools) 
devtools::install_github("aplantin/pldist")

## ----load----------------------------------------------------------------
library(pldist)

## ----load-data-----------------------------------------------------------
data("sim.tree")
data("paired.otus"); data("paired.meta")
data("bal.long.otus"); data("bal.long.meta")
data("unbal.long.otus"); data("unbal.long.meta")

## ----transform-----------------------------------------------------------
# Input: Notice that row names are sample IDs 
paired.otus[1:4,1:4]
paired.meta[1:4,]

# Transformation function 
res <- pltransform(paired.otus, paired.meta, paired = TRUE, check.input = TRUE)

# Binary transformation 
# 0.5 indicates OTU was present at Time 2, absent at Time 1
# -0.5 indicates OTU was present at Time 1, absent at Time 2 
# Row names are now subject IDs 
res$dat.binary   

# Quantitative transformation (see details in later sections)
round(res$dat.quant, 2)

# Average proportion per OTU per subject 
round(res$avg.prop, 2)

# This was a paired transformation 
res$type 

# For comparison, this uses a longitudinal transformation (applied at 2 time points)
# due to the argument "paired = FALSE". 
# Type is "Balanced" because same time points were observed for all subjects 
res2 <- pltransform(paired.otus, paired.meta, paired = FALSE, check.input = TRUE)
res2$type   

# With the longitudinal binary transformation applied at 2 time points, the value 
# is 1 if any change in presence/absence was observed, 0 otherwise 
res2$dat.binary 

# And if you use an unbalanced design, the function gives a warning. 
# It otherwise operates in the same way. 
res3 <- pltransform(unbal.long.otus, unbal.long.meta, paired = FALSE, check.input = TRUE)
round(res3$dat.quant[1:4,1:4], 2)


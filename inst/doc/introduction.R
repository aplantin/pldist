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

## ----lunifrac------------------------------------------------------------
# Input: 
#    otu.tab: OTU matrix (as above) 
#    metadata: Metadata (as above) 
#    tree: Rooted phylogenetic tree with tip labels that match OTU column names 
#    gam: Gamma parameters, vector of values between 0 and 1
#         (controls weight placed on abundant OTUs) 
#    paired: Indicates whether to use paired transformation (TRUE) or longitudinal (FALSE) 
#    check.input: Indicates whether to check function input before proceeding (default TRUE)
#    
D.unifrac <- LUniFrac(otu.tab = paired.otus, metadata = paired.meta, tree = sim.tree, 
                      gam = c(0, 0.5, 1), paired = TRUE, check.input = TRUE)
D.unifrac[, , "d_1"]   # gamma = 1 (quantitative paired transformation)
D.unifrac[, , "d_UW"]  # unweighted LUniFrac (qualitative/binary paired transf.)


# Same procedure for longitudinal data 
D2.unifrac <- LUniFrac(otu.tab = bal.long.otus, metadata = bal.long.meta, tree = sim.tree, 
                      gam = c(0, 0.5, 1), paired = FALSE, check.input = TRUE)
D2.unifrac[, , "d_1"]   # gamma = 1 (quantitative longitudinal transformation)
D2.unifrac[, , "d_UW"]  # unweighted LUniFrac (qualitative/binary longitudinal transf.)

## ----pldist--------------------------------------------------------------
# Gower distance, paired & quantitative transformation 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "gower")$D

# Gower distance, paired & qualitative/binary transformation 
pldist(paired.otus, paired.meta, paired = TRUE, binary = TRUE, method = "gower")$D

# Gower distance, longitudinal & quantitative transformation 
pldist(bal.long.otus, bal.long.meta, paired = FALSE, binary = FALSE, method = "gower")$D

# Gower distance, longitudinal & qualitative/binary transformation 
pldist(bal.long.otus, bal.long.meta, paired = FALSE, binary = TRUE, method = "gower")$D

# Other distances 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "bray")$D
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "kulczynski")$D
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "jaccard")$D

# UniFrac additionally requires a phylogenetic tree and gamma values 
# (Gamma controls weight placed on abundant lineages) 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, 
    method = "unifrac", tree = sim.tree, gam = c(0, 0.5, 1))$D 

## ----gen-tree------------------------------------------------------------
# tree tip names must match column names in OTU table
gen.tree <- function(seed, notus) {
  set.seed(seed)
  sim.tree = rtree(n=notus)
  sim.tree$tip.label <- paste("otu", sample(1:notus), sep = "")
  return(sim.tree)
}

## ----data-gen-fxn--------------------------------------------------------
# Parameters: 
#     nsubj: number of subjects 
#     maxtimes: maximum number of time points 
#         (use maxtimes=2 for paired data) 
#     maxdiff: maximum difference between observed times 
#         (use maxdiff = 1 for observation at consecutive time units) 
#     balanced: logical indicating whether design should be balanced 
#         (all subjects observed at the same time points)
#     notus: number of observed OTUs 
#     propzero: proportion of table cells that should be zero
#         (microbiome data tends to have a high proportion of zeros) 
#     maxct: maximum possible read count in a single cell
gen.data <- function(seed, nsubj, maxtimes, maxdiff, balanced, notus, propzero, maxct) {
  set.seed(seed)
  if (maxtimes == 2) {
    ntimes <- rep(2, nsubj) 
  } else {
    if (balanced) {
      ntimes <- rep(sample(2:maxtimes)[1], nsubj)
    } else {
      ntimes <- sample(2:maxtimes, nsubj, replace = TRUE)
    }
  }
  ncells = sum(ntimes) * notus
  nzero = floor(ncells*propzero)
  toy.otus <- matrix(0, nrow = sum(ntimes), ncol = notus)
  while (any(c(apply(toy.otus, 1, FUN = function(x) all(x == 0)), 
               apply(toy.otus, 2, FUN = function(x) all(x == 0))))) {
    toy.otus <- matrix(sample(c(sample(1:maxct, (ncells - nzero), replace = TRUE), rep(0, nzero))), 
                       nrow = sum(ntimes), ncol = notus)
  }
  toy.props <- counts2props(toy.otus) 
  subjIDs <- unlist(sapply(1:nsubj, FUN = function(i) {
    rep(paste("subj", i, sep = ""), ntimes[i]) }, simplify = FALSE))
  sampIDs <- paste(unlist(sapply(1:nsubj, FUN = function(i) {
    rep(paste("subj", i, sep = ""), ntimes[i]) }, simplify = FALSE)), 
    unlist(sapply(1:nsubj, FUN = function(i) letters[1:ntimes[i]], simplify = FALSE)), sep = "")
  if (balanced) {
    times <- rep(cumsum(c(1, sample(1:maxdiff, ntimes[1]-1, replace = TRUE))), nsubj)
  } else {
    times <- unlist(sapply(1:nsubj, FUN = function(i) {
      cumsum(c(1, sample(1:maxdiff, ntimes[i]-1, replace = TRUE)))}, simplify = FALSE))
  }
  toy.meta <- data.frame(subjID = subjIDs, sampID = sampIDs, 
                         time   = times, stringsAsFactors = FALSE)
  rownames(toy.otus) = toy.meta$sampID
  colnames(toy.otus) = paste("otu", 1:notus, sep = "")
  return(list(otus = toy.otus, metadata = toy.meta))
}

## ----gen-test-data, eval=FALSE-------------------------------------------
#  # Paired data (two sequential observations on each subject)
#  paired.data <- gen.data(seed = 1, nsubj = 5, maxtimes = 2, maxdiff = 1, balanced = TRUE, notus = 10, propzero = 0.5, maxct = 500)
#  paired.otus <- paired.data$otus
#  paired.meta <- paired.data$metadata
#  
#  # Balanced longitudinal data: Same number of observations, at same times, for each subject.
#  # Here each subject has three observations total, at days 1, 3, and 6.
#  # (up to 4 observations allowed, depending on random seed, with max 3 time units between observations)
#  bal.long.data <- gen.data(seed = 4, nsubj = 5, maxtimes = 4, maxdiff = 3, balanced = TRUE, notus = 10, propzero = 0.5, maxct = 500)
#  bal.long.otus <- bal.long.data$otus
#  bal.long.meta <- bal.long.data$metadata
#  
#  # Unbalanced longitudinal data: Different number and spacing of observations for subjects.
#  # (up to 4 observations per subject with up to 3 time units between observations)
#  unbal.long.data <- gen.data(seed = 1, nsubj = 5, maxtimes = 4, maxdiff = 3, balanced = FALSE, notus = 10, propzero = 0.5, maxct = 500)
#  unbal.long.otus <- unbal.long.data$otus
#  unbal.long.meta <- unbal.long.data$metadata
#  
#  # Rooted phylogenetic tree with 10 OTUs
#  sim.tree <- gen.tree(seed = 1, notus = 10)


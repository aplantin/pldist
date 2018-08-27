#' Simulated OTU data for paired study design. 
#'
#' Simulation code is included in the package vignette. 
#' Corresponding metadata is stored in `paired.meta`. 
#' 
#' @usage data(paired.otus)
#'
#' @format A matrix with 10 rows and 10 columns. Rows are samples, columns are OTUs. 
"paired.otus"


#' Simulated metadata for paired study design. 
#'
#' Simulation code is included in the package vignette. 
#' Corresponding OTU matrix is stored in `paired.otus`. 
#' 
#' @usage data(paired.meta)
#'
#' @format A data frame with 10 rows and 3 columns. 
#' \describe{
#'   \item{subjID}{Subject identifiers}
#'   \item{sampID}{Sample identifiers, matches row names of OTU count matrix}
#'   \item{time}{Time indicator, takes values 1 or 2}
#' }
"paired.meta"


#' Simulated OTU data for balanced longitudinal study design. 
#'
#' Simulation code is included in the package vignette. 
#' Corresponding metadata is stored in `bal.long.meta`. 
#' 
#' @usage data(bal.long.otus)
#'
#' @format A matrix with 15 rows and 10 columns. Rows are samples, columns are OTUs. 
"bal.long.otus"


#' Simulated metadata for balanced longitudinal study design. 
#'
#' Simulation code is included in the package vignette. 
#' Corresponding OTU matrix is stored in `bal.long.otus`. 
#' 
#' @usage data(bal.long.meta)
#'
#' @format A data frame with 15 rows and 3 columns. 
#' \describe{
#'   \item{subjID}{Subject identifiers}
#'   \item{sampID}{Sample identifiers, matches row names of OTU count matrix}
#'   \item{time}{Time indicator}
#' }
"bal.long.meta"

#' Simulated OTU data for unbalanced longitudinal study design. 
#'
#' Simulation code is included in the package vignette. 
#' Corresponding metadata is stored in `unbal.long.meta`. 
#' 
#' @usage data(unbal.long.otus)
#'
#' @format A matrix with 14 rows and 10 columns. Rows are samples, columns are OTUs. 
"unbal.long.otus"


#' Simulated metadata for balanced longitudinal study design. 
#'
#' Simulation code is included in the package vignette. 
#' Corresponding OTU matrix is stored in `unbal.long.otus`. 
#'
#' @format A data frame with 14 rows and 3 columns. 
#' \describe{
#'   \item{subjID}{Subject identifiers}
#'   \item{sampID}{Sample identifiers, matches row names of OTU count matrix}
#'   \item{time}{Time indicator}
#' }
#' 
#' @usage data(unbal.long.meta)
#' 
"unbal.long.meta"


#' Simulated rooted phylogenetic tree. 
#'
#' Simulation code is included in the package vignette. 
#' Tree includes 10 OTUs and may be used with any of the 
#' simulated data sets (paired, balanced longitudinal, or 
#' unbalanced longitudinal). 
#' 
#' @usage data(sim.tree)
#'
#' @format An object of class "phylo". 
#' 
"sim.tree"


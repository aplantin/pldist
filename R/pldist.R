#' pldist 
#'
#' Function that calculates paired and longitudinal ecological distance/dissimilarity 
#' matrices. Includes qualitative and quantitative versions of Bray-Curtis, Jaccard, Kulczynski, 
#' Gower, and unweighted and generalized UniFrac distances/dissimilarities. UniFrac-based 
#' metrics are based in part on GUniFrac (Jun Chen & Hongzhe Li (2012)).  
#'
#' @param otus OTU count or frequency table, containing one row per sample and one column per OTU. 
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier (if using longitudinal distances, this must be numeric or 
#'     convertable to numeric).
#' @param paired Logical indicating whether to use the paired version of the metric (TRUE) or the 
#'     longitudinal version (FALSE). Paired analyis is only possible when there are exactly 2 
#'     unique time points/identifiers for each subject or pair. 
#' @param binary Logical indicating whether to use the qualitative (TRUE) or quantitative (FALSE) 
#'     version of each metric. Qualitative analysis only incorporates changes in OTU presence or 
#'     absence; quantitative analysis incorporates changes in abundance. 
#' @param method Desired distance metric. Choices are braycurtis, jaccard, kulczynski, gower, and 
#'     unifrac, or any unambiguous abbreviation thereof. 
#' @param tree Rooted phylogenetic tree of R class "phylo". Default NULL; only needed for 
#'     UniFrac family distances. 
#' @param gam Parameter controlling weight on abundant lineages for UniFrac family distances. The 
#'     same weight is used within a subject as between subjects. Default (0, 0.5, 1). 
#' @return Returns a list with elements: 
#'     \item{D}{If any metric other than UniFrac is used, D is an n x n distance (or dissimilarity) matrix. 
#'     For UniFrac-family dissimilarities, D is a (K+1) dimensional array containing the paired or 
#'     longitudinal UniFrac dissimilarities with the K specified gamma values plus the unweighted 
#'     distance. The unweighted distance matrix may be accessed by result[,,"d_UW"], and the 
#'     generalized dissimilarities by result[,,"d_G"] where G is the particular choice of gamma.} 
#'     \item{type}{String indicating what type of dissimilarity was requested.}
#'
#' @export
#'
pldist <- function(otus, metadata, paired = FALSE, binary = FALSE, method, tree = NULL, gam = c(0, 0.5, 1)) {
  ## Find desired method 
  method.opts = c("braycurtis", "jaccard", "kulczynski", "gower", "unifrac")
  this.method = pmatch(trimws(tolower(method)), method.opts, nomatch = NA)
  if (is.na(this.method)) {
    stop("Method does not match any expected methods. Please see list of options in documentation.")
  } 
  method = method.opts[this.method]
  
  ## Check and prepare all input 
  if (nrow(otus) != nrow(metadata)) {
    stop("Number of rows of metadata and OTUs should match") } 
  
  if (any(apply(otus, 1, FUN = function(x) all(x == 0)))) {
    stop("At least one subject has uniformly zero OTU counts. Please exclude.") }
  
  if (any(apply(otus, 2, FUN = function(x) all(x == 0)))) {
    warning("Some OTUs have count zero for all subjects and are being excluded.")
    otus <- otus[, which(apply(otus, 2, FUN = function(x) !all(x == 0)))]
  }
  
  if (!all(colnames(metadata) == c("subjID", "sampID", "time"))) {
    stop("Please format metadata with columns \"subjID\", \"sampID\", \"time\" in that order") }
  
  if (is.null(rownames(otus)) | !all(rownames(otus) == metadata$sampID)) {
    stop("Please ensure rownames of OTU matrix exactly match sample IDs in metadata") }
  
  if (!all(apply(otus, 1, sum) == 1))  otus <- counts2props(otus) 
  
  if (paired) { 
    if (length(unique(metadata$time)) > 2) {
      stop("Paired dissimilarities were requested, but >2 unique time points/groups were provided.")
    } else if (length(unique(metadata$time)) < 2) {
      stop("Paired dissimilarities were requested, but <2 unique time points/groups were provided.")
    } 
    
    persubj <- aggregate(metadata$time, by = list(metadata$subjID), FUN = function(x) length(x))$x
    if (length(unique(persubj)) != 1) {
      stop("Paired dissimilarities were requested, but some groups/subjects do not have 2 observations. \n
           Please check for missing or miscoded data and exclude any unpaired observations.")
    }
    metadata$time = as.numeric(as.factor(metadata$time))
  } else {      
    metadata$time = as.numeric(metadata$time) 
  }
  
  ## Calculate distances/dissimilarities 
  
  if (method.opts[this.method] != "unifrac") {
    ## Calculate transformed data and apply distance (all except UniFrac) 
    tsf.res <- pl.transform(otus = otus, metadata = metadata, paired = paired)
    D <- switch(method, 
                braycurtis = braycurtis(tsf.res$tsf.data, binary = binary), 
                jaccard = jaccard(tsf.res$tsf.data, paired = paired, binary = binary), 
                kulczynski = kulczynski(tsf.res$tsf.data, paired = paired, binary = binary), 
                gower = gower(tsf.res$tsf.data, binary = binary)
                ) 
  } else {
    ## Calculate paired/longitudinal UniFrac dissimilarities 
    if (is.null(tree)) stop("Tree is required for UniFrac family metrics.")
    if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!") 
    if (paired) {
      D <- PUniFrac(otu.tab = otus, tree = tree, gam = gam, metadata = metadata)
    } else {
      D <- LUniFrac(otu.tab = otus, tree = tree, gam = gam, metadata = metadata)
    }
  } 
  
  if (paired) {
    if (method != "unifrac") {
      if (binary) {
        type = paste("Method: ", method, "; Paired, Binary") 
      } else {
        type = paste("Method: ", method, "; Paired, Quantitative") 
      }
    } else {
      type = paste("UniFrac family; Paired") 
    } 
  } else {
    if (method != "unifrac") {
      if (binary) {
        type = paste("Method: ", method, "; Longitudinal, Binary") 
      } else {
        type = paste("Method: ", method, "; Longitudinal, Quantitative") 
      }
    } else {
      type = paste("UniFrac family; Longitudinal") 
    } 
  }
  return(list(D = D, type = type)) 
}
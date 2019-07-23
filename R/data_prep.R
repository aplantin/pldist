#' data_prep
#' 
#' Checks data input for transformation (pltransform) and dissimilarity (pldist) functions. Prepares 
#' matrix of OTU proportions and CLR-transformed matrix. 
#' 
#' @param otus Matrix of OTU counts (better) or proportions. Row names must be sample identifiers 
#'     (matching metadata), and column names must be OTU identifiers (enforced if 
#'     using UniFrac distances). 
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier (if using longitudinal distances, this must be numeric or 
#'     convertable to numeric). 
#' @param paired Logical indicating whether to use the paired version of the metric (TRUE) or the 
#'     longitudinal version (FALSE). Paired analyis is only possible when there are exactly 2 
#'     unique time points/identifiers for each subject or pair. 
#' @param pseudoct Pseudocount value to be added to each cell of the matrix. Default is NULL; if NULL, 
#'     0.5 will be added if data are counts, min(1e-06, 0.5*min(nonzero p)) will be added if data are 
#'     proportions, and nothing will be added if no cells have zero values. This is only done for the 
#'     CLR-transformed data matrix, not the OTU proportion matrix. 
#' @return Returns OTU proportions, CLR-transformed OTU proportions, and metadata files, all checked 
#'     for formatting and value problems. 
#'
#' @export
#' 
data_prep <- function(otus, metadata, paired, pseudoct = NULL) {
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
  
  if (!all(apply(otus, 1, sum) == 1))  {
    otu.props <- counts2props(otus) 
  } else {
    otu.props <- otus 
  }
  otu.clr <-  psct_clr(otus, pseudocount = pseudoct)
  rownames(otu.props) = rownames(otu.clr) = rownames(otus)
  colnames(otu.props) = colnames(otu.clr) = colnames(otus)
  
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
    persubj <- aggregate(metadata$time, by = list(metadata$subjID), FUN = function(x) length(x))$x
    if (any(persubj < 2)) {
      stop("Some group(s) or subject(s) do not have at least 2 observations. \n
           Please check for missing or miscoded data and exclude singleton groups/subjects.")
    }
    metadata$time = as.numeric(metadata$time) 
  }
  
  return(list(otu.props = otu.props, otu.clr = otu.clr, metadata = metadata))
}

#' psct_clr
#' 
#' Adds a pseudocount, (re-)closes to proportions, and then applies centered log-ratio transformation 
#' 
#' @param otus Matrix of OTU counts or proportions. Rows are samples, columns are OTUs. 
#' @param pseudocount Pseudocount to be added to all values in OTU matrix. If NULL, then 
#'     0.5 is added if entries are counts, or min(1e-06, 0.5*min(nonzero p)) if entries
#'     are proportions. If all entries are nonzero, nothing is added. 
#' @return Matrix of CLR-transformed data (dimensions, rownames, colnames match input)
#' 
#' @export 
psct_clr <- function(otus, pseudocount = NULL) {
  if (any(otus == 0)) {
    ## pseudocount 
    if (is.null(pseudocount)) {
      subjsums <- apply(otus, 1, FUN = function(x) sum(x))
      if (all(subjsums == 1)) {  # data are proportions 
        otus <- otus + min(1e-06, otus[otus != 0])
      } else {
        otus <- otus + 0.5
      }
    } else {
      otus <- otus + pseudocount 
    }
  }
  
  otus <- counts2props(otus) 
  
  gm <- apply(otus, 1, FUN = function(x) exp(mean(log(x))))
  tsf.otus <- matrix(nrow = nrow(otus), ncol = ncol(otus)) 
  for (i in 1:nrow(otus)) {
    tsf.otus[i,] = log(otus[i,]/gm[i])
  }
  rownames(tsf.otus) = rownames(otus) 
  colnames(tsf.otus) = colnames(otus) 
  return(tsf.otus)
}


#' counts2props 
#' 
#' Converts OTU counts to OTU proportions/relative abundances. 
#'     
#' @param x Matrix of OTU counts (rows are subjects, columns are taxa). 
#' 
#' @return n x p matrix of OTU proportions.  
#' 
#' @export 
#' 
counts2props <- function(x) {
  mat <- t(apply(x, 1, FUN = function(y) y/sum(y)))
  rownames(mat) = rownames(x)
  colnames(mat) = colnames(x)
  return(mat)
}



check_input <- function(otus, metadata, paired) {
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
  
  return(list(otus = otus, metadata = metadata))
}
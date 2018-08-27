#' tsf_paired 
#' 
#' OTU transformation for paired data. Computes within-subject change (in presence 
#'     for qualitative metrics and abundance for quantitative metrics) between time 
#'     points for each taxon. 
#'     
#' @param otus Matrix of OTU counts or proportions. Notes: (1) Will be transformed to 
#'     proportions if it's not already; (2) Row names must be sample identifiers 
#'     (matching metadata), and column names must be OTU identifiers (enforced if 
#'     using UniFrac distances). 
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier (must have two unique values for paired transformation). 
#' @importFrom stats aggregate
#' 
#' @return List with the following elements. Both data matrices have subject identifiers 
#'     as row names and OTU identifiers as column names.  
#'     \item{dat.binary}{n x p matrix of data after paired, binary/qualitative transformation} 
#'     \item{dat.quant}{n x p matrix of data after paired, quantitative transformation} 
#'     \item{avg.prop}{n x p matrix with overall average proportion of each taxon} 
#' 
tsf_paired <- function(otus, metadata) {
  ## Prepare output data frame 
  n <- length(unique(metadata$subjID))
  out.data <- matrix(0, nrow = n, ncol = ncol(otus))
  rownames(out.data) <- unique(metadata$subjID)
  colnames(out.data) <- colnames(otus) 
  
  ## Main function 
  out.binary = out.quant = out.avgprop = out.data 
  for (i in 1:nrow(out.data)) {
    t1.idx <- which(metadata$subjID == rownames(out.data)[i] & metadata$time == 1)
    t2.idx <- which(metadata$subjID == rownames(out.data)[i] & metadata$time == 2)
    out.binary[i, ] <- 0.5 * (as.numeric(otus[t2.idx,] > 0) - as.numeric(otus[t1.idx,] > 0)) 
    nonz <- which(otus[t2.idx,] != 0 | otus[t1.idx,] != 0) 
    out.quant[i, nonz] <- 0.5 * (otus[t2.idx, nonz] - otus[t1.idx, nonz]) / (otus[t2.idx, nonz] + otus[t1.idx, nonz]) 
    out.avgprop[i, ] <- 0.5 * (otus[t2.idx, ] + otus[t1.idx, ])
  }
  return(list(dat.binary = out.binary, dat.quant = out.quant, avg.prop = out.avgprop))   
} 

#' tsf_paired 
#' 
#' OTU transformation for longitudinal data. Computes average within-subject change 
#'     (in presence for qualitative metrics, abundance for quantitative metrics) 
#'     during one unit of time for each taxon. 
#'     
#' @param otus Matrix of OTU counts or proportions. Notes: (1) Will be transformed to 
#'     proportions if it's not already; (2) Row names must be sample identifiers 
#'     (matching metadata), and column names must be OTU identifiers (enforced if 
#'     using UniFrac distances). 
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier (if using longitudinal distances, this must be numeric or 
#'     convertable to numeric). 
#' 
#' @return List with the following elements. Both data matrices have subject identifiers 
#'     as row names and OTU identifiers as column names.  
#'     \item{dat.binary}{n x p matrix of data after longitudinal, binary/qualitative transformation} 
#'     \item{dat.quant}{n x p matrix of data after longitudinal, quantitative transformation} 
#'     \item{avg.prop}{n x p matrix with overall average proportion of each taxon} 
#' 
tsf_long <- function(otus, metadata) {
  ## Prepare output data frame
  n <- length(unique(metadata$subjID))
  out.data <- matrix(0, nrow = n, ncol = ncol(otus))
  rownames(out.data) <- unique(metadata$subjID)
  colnames(out.data) <- colnames(otus) 
  
  ## Main function 
  out.binary = out.quant = out.data 
  out.avgprop = out.data 
  
  for (i in 1:nrow(out.data)) {
    ## Prep subject 
    subj.idx <- which(metadata$subjID == rownames(out.data)[i])
    subj.otu <- otus[subj.idx, ]
    subj.times <- metadata$time[subj.idx] 
    
    ord <- order(metadata$time[subj.idx])
    subj.otu <- subj.otu[ord, ]
    subj.times <- subj.times[ord]
    qi <- nrow(subj.otu)
    
    ## Calculate both 
    dk.uw <- rep(0, ncol(otus))
    dk.g <- rep(0, ncol(otus))
    cumprop <- subj.otu[1,] 
    for (j in 1:(qi-1)) {
      dk.uw = dk.uw + (1/(subj.times[j+1] - subj.times[j])) * abs(as.numeric(subj.otu[(j+1), ] > 0) - as.numeric(subj.otu[j, ] > 0))
      nonz <- which(subj.otu[(j+1), ] != 0 | subj.otu[j, ] != 0)
      dk.g[nonz] = dk.g[nonz] + (1/(subj.times[j+1] - subj.times[j])) * 
        abs((subj.otu[(j+1), nonz] - subj.otu[j, nonz])/(subj.otu[(j+1), nonz] + subj.otu[j, nonz]))
      cumprop = cumprop + subj.otu[(j+1), ]
    }
    dk.uw = dk.uw/(qi - 1)
    dk.g = dk.g/(qi - 1)
    cumprop = cumprop/qi 
    
    ## Fill row 
    out.binary[i, ] <- dk.uw 
    out.quant[i, ] <- dk.g 
    out.avgprop[i, ] <- cumprop 
  }
  return(list(dat.binary = out.binary, dat.quant = out.quant, avg.prop = out.avgprop)) 
}

#' counts2props 
#' 
#' Converts OTU counts to OTU proportions/relative abundances. 
#'     
#' @param x Matrix of OTU counts (rows are subjects, columns are taxa). 
#' 
#' @return n x p matrix of OTU proportions.  
#' 
counts2props <- function(x) {
  return(t(apply(x, 1, FUN = function(y) y/sum(y))))
}





#' tsf_paired 
#' 
#' OTU transformation for longitudinal data. Computes average within-subject change 
#'     (in presence for qualitative metrics, abundance for quantitative metrics) 
#'     during one unit of time for each taxon. 
#'     
#' @param otus Matrix of OTU counts or proportions. Notes: (1) Will be transformed to 
#'     proportions if it's not already; (2) Row names must be sample identifiers 
#'     (matching metadata), and column names must be OTU identifiers (enforced if 
#'     using UniFrac distances). 
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier (if using longitudinal distances, this must be numeric or 
#'     convertable to numeric). 
#' @param paired Logical indicating whether to use the paired version of the metric (TRUE) or the 
#'     longitudinal version (FALSE). Paired analyis is only possible when there are exactly 2 
#'     unique time points/identifiers for each subject or pair. 
#' @param check.input Logical indicating whether to check input values (default TRUE). 
#' 
#' @return List with the following elements. Both data matrices have subject identifiers 
#'     as row names and OTU identifiers as column names.  
#'     \item{tsf.data}{List with 3 elements: 
#'         (1) dat.binary: n x p matrix of data after longitudinal, binary/qualitative transformation 
#'         (2) dat.quant: n x p matrix of data after longitudinal, quantitative transformation
#'         (3) avg.prop: n x p matrix with overall average proportion of each taxon }
#'     \item{type}{Type of transformation that was used (paired, balanced longitudinal, 
#'     unbalanced longitudinal) with a warning if unbalanced longitudinal.} 
#'     
#' 
pltransform <- function(otus, metadata, paired, check.input = TRUE) {
  if (check.input) {
    okdat <- check_input(otus, metadata, paired)
    otus <- okdat$otus 
    metadata <- okdat$metadata 
    remove(okdat) 
  }
  
  ## calculate appropriate transformations 
  if (paired) {
    res <- tsf_paired(otus, metadata)
  } else {
    res <- tsf_long(otus, metadata)
    if (length(unique(table(metadata$time))) != 1) { balanced = FALSE } else { balanced = TRUE } 
  }
  if (paired) { type = "paired" 
  } else if (balanced) { type = "balanced longitudinal"
  } else { 
    type = "unbalanced longitudinal (WARNING: this transformation is not recommended for strongly unbalanced designs!)" 
    warning("WARNING: this transformation is not recommended for strongly unbalanced designs!")
    }
   
  ## return 
  return(list(tsf.data = res, type = type))
}



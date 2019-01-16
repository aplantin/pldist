#' tsf_paired 
#' 
#' OTU transformation for paired data. Computes within-subject change (in presence 
#'     for qualitative metrics and abundance for quantitative metrics) between time 
#'     points for each taxon. 
#'     
#' @param otu.props Matrix of OTU proportions (rownames = sample IDs, colnames = OTU IDs) 
#' @param otu.clr Matrix of CLR-transformed OTU proportions (rownames = sample IDs, colnames = OTU IDs)
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier (must have two unique values for paired transformation). 
#' @param norm Indicator of whether to normalize the difference to average taxon abundance or not (default TRUE)
#' @importFrom stats aggregate
#' 
#' @return List with the following elements. Data matrices have subject identifiers 
#'     as row names and OTU identifiers as column names.  
#'     \item{dat.binary}{n x p matrix of data after paired, binary/qualitative transformation} 
#'     \item{dat.quant.prop}{n x p matrix of data after paired, quantitative transformation applied to OTU proportions} 
#'     \item{dat.quant.clr}{n x p matrix of data after paired, quantitative transformation applied to CLR-transformed OTU proportions} 
#'     \item{avg.prop}{n x p matrix with overall average proportion of each taxon}
#'     
#' @export
#' 
tsf_paired <- function(otu.props, otu.clr, metadata, norm = TRUE) {
  ## Prepare output data frame
  n <- length(unique(metadata$subjID))
  out.data <- matrix(0, nrow = n, ncol = ncol(otu.props))
  rownames(out.data) <- unique(metadata$subjID)
  colnames(out.data) <- colnames(otu.props)
  
  ## Main function
  out.binary = out.quant.prop = out.quant.clr = out.avgprop = out.data
  for (i in 1:nrow(out.data)) {
    t1.idx <- which(metadata$subjID == rownames(out.data)[i] & metadata$time == 1)
    t2.idx <- which(metadata$subjID == rownames(out.data)[i] & metadata$time == 2)
    out.binary[i, ] <- 0.5 * (as.numeric(otu.props[t2.idx,] > 0) - as.numeric(otu.props[t1.idx,] > 0))
    nonz <- which(otu.props[t2.idx,] != 0 | otu.props[t1.idx,] != 0)
    nonz2 <- which(otu.clr[t2.idx,] != 0 | otu.clr[t1.idx,] != 0)
    if (norm) {
      out.quant.prop[i, nonz] <- 0.5 * (otu.props[t2.idx, nonz] - otu.props[t1.idx, nonz]) / (otu.props[t2.idx, nonz] + otu.props[t1.idx, nonz])
      out.quant.clr[i, nonz2] <- 0.5 * (otu.clr[t2.idx, nonz2] - otu.clr[t1.idx, nonz2]) / (abs(otu.clr[t2.idx, nonz2]) + abs(otu.clr[t1.idx, nonz2]))
    } else {
      out.quant.prop[i, nonz] <- 0.5 * (otu.props[t2.idx, nonz] - otu.props[t1.idx, nonz])
      out.quant.clr[i, nonz2] <- 0.5 * (otu.clr[t2.idx, nonz2] - otu.clr[t1.idx, nonz2])
    }
    out.avgprop[i, ] <- 0.5 * (otu.props[t2.idx, ] + otu.props[t1.idx, ])
  }
  return(list(dat.binary = out.binary, dat.quant.prop = out.quant.prop, dat.quant.clr = out.quant.clr, avg.prop = out.avgprop))
}


#' tsf_long 
#' 
#' OTU transformation for longitudinal data. Computes average within-subject change 
#'     (in presence for qualitative metrics, abundance for quantitative metrics) 
#'     during one unit of time for each taxon. 
#'     
#' @param otu.props Matrix of OTU proportions (rownames = sample IDs, colnames = OTU IDs) 
#' @param otu.clr Matrix of CLR-transformed OTU proportions (rownames = sample IDs, colnames = OTU IDs)
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier.
#' @param norm Indicator of whether to normalize the difference to average taxon abundance or not (default TRUE)
#' 
#' @return List with the following elements. Data matrices have subject identifiers 
#'     as row names and OTU identifiers as column names.  
#'     \item{dat.binary}{n x p matrix of data after longitudinal, binary/qualitative transformation} 
#'     \item{dat.quant.prop}{n x p matrix of data after longitudinal, quantitative transformation applied to OTU proportions} 
#'     \item{dat.quant.clr}{n x p matrix of data after longitudinal, quantitative transformation applied to CLR-transformed OTU proportions} 
#'     \item{avg.prop}{n x p matrix with overall average proportion of each taxon}
#'     
#' @export
#'  
tsf_long <- function(otu.props, otu.clr, metadata, norm = TRUE) {
  ## Prepare output data frame
  n <- length(unique(metadata$subjID))
  out.data <- matrix(0, nrow = n, ncol = ncol(otu.props))
  rownames(out.data) <- unique(metadata$subjID)
  colnames(out.data) <- colnames(otu.props)
  
  ## Main function
  out.binary = out.quant.props = out.quant.clr = out.avgprop = out.data
  
  for (i in 1:nrow(out.data)) {
    ## Prep subject
    subj.idx <- which(metadata$subjID == rownames(out.data)[i])
    subj.otu.props <- otu.props[subj.idx, ]
    subj.otu.clr <- otu.clr[subj.idx, ]
    subj.times <- metadata$time[subj.idx]
    
    ord <- order(metadata$time[subj.idx])
    subj.otu.props <- subj.otu.props[ord, ]
    subj.otu.clr <- subj.otu.clr[ord, ]
    subj.times <- subj.times[ord]
    qi <- nrow(subj.otu.props)
    
    ## Calculate all
    dk.uw <- rep(0, ncol(otu.props))
    dk.gp <- rep(0, ncol(otu.props))
    dk.gclr <- rep(0, ncol(otu.clr))
    cumprop <- subj.otu.props[1,]
    for (j in 1:(qi-1)) {
      dk.uw = dk.uw + (1/(subj.times[j+1] - subj.times[j])) *
        abs(as.numeric(subj.otu.props[(j+1), ] > 0) - as.numeric(subj.otu.props[j, ] > 0))
      nonz <- which(subj.otu.props[(j+1), ] != 0 | subj.otu.props[j, ] != 0)
      nonz2 <- which(subj.otu.clr[(j+1), ] != 0 | subj.otu.clr[j, ] != 0)
      if (norm) {
        dk.gp[nonz] = dk.gp[nonz] + (1/(subj.times[j+1] - subj.times[j])) *
          abs((subj.otu.props[(j+1), nonz] - subj.otu.props[j, nonz])/(subj.otu.props[(j+1), nonz] + subj.otu.props[j, nonz]))
        dk.gclr[nonz2] = dk.gclr[nonz2] + (1/(subj.times[j+1] - subj.times[j])) *
          abs(subj.otu.clr[(j+1), nonz2] - subj.otu.clr[j, nonz2])/(abs(subj.otu.clr[(j+1), nonz2]) + abs(subj.otu.clr[j, nonz2]))
      } else {
        dk.gp[nonz] = dk.gp[nonz] + (1/(subj.times[j+1] - subj.times[j])) *
          abs(subj.otu.props[(j+1), nonz] - subj.otu.props[j, nonz])
        dk.gclr[nonz2] = dk.gclr[nonz2] + (1/(subj.times[j+1] - subj.times[j])) *
          abs(subj.otu.clr[(j+1), nonz2] - subj.otu.clr[j, nonz2])
      }
      cumprop = cumprop + subj.otu.props[(j+1), ]
    }
    dk.uw = dk.uw/(qi - 1)
    dk.gp = dk.gp/(qi - 1)
    dk.gclr = dk.gclr/(qi - 1)
    cumprop = cumprop/qi
    
    ## Fill row
    out.binary[i, ] <- dk.uw
    out.quant.props[i, ] <- dk.gp
    out.quant.clr[i, ] <- dk.gclr
    out.avgprop[i, ] <- cumprop
  }
  return(list(dat.binary = out.binary, dat.quant.props = out.quant.props, dat.quant.clr = out.quant.clr, avg.prop = out.avgprop))
}


#' pltransform 
#' 
#' OTU transformation for paired and longitudinal data. Computes average within-subject 
#'     change (in presence for qualitative metrics, proportional or CLR-transformed 
#'     abundance for quantitative metrics) during one unit of time for each taxon. 
#'     
#' @param otu.data OTU data after pre-processing using data_prep() function. List with elements 
#'     otu.props, otu.clr, and metadata. 
#' @param paired Logical indicating whether to use the paired version of the metric (TRUE) or the 
#'     longitudinal version (FALSE). Paired analyis is only possible when there are exactly 2 
#'     unique time points/identifiers for each subject or pair. 
#' @param norm Logical indicating whether quantitative differences should be normalized by 
#'     taxon abundance (default TRUE)
#' 
#' @return List with the following elements. Data matrices have subject identifiers 
#'     as row names and OTU identifiers as column names.   
#'     \item{dat.binary}{n x p matrix of data after binary/qualitative transformation} 
#'     \item{dat.quant.prop}{n x p matrix of data after quantitative transformation applied to proportions}
#'     \item{dat.quant.prop}{n x p matrix of data after quantitative transformation applied to CLR-transformed proportions}
#'     \item{avg.prop}{n x p matrix with overall average proportion of each taxon}
#'     \item{type}{Type of transformation that was used (paired, balanced longitudinal, unbalanced longitudinal) with a warning if unbalanced longitudinal.} 
#' @examples
#' data("paired.otus")
#' data("paired.meta")
#'  # paired transformation
#' otudat <- data_prep(paired.otus, paired.meta, paired = TRUE) 
#' res1 <- pltransform(otudat, paired = TRUE, norm = TRUE) 
#'  # longitudinal transformation 
#' otudat <- data_prep(paired.otus, paired.meta, paired = FALSE)
#' res2 <- pltransform(otudat, paired = FALSE, norm = TRUE) 
#'     
#' @export 
#' 
pltransform <- function(otu.data, paired, norm = TRUE) {
  if (!all(names(otu.data) == c("otu.props", "otu.clr", "metadata"))) {
    stop("Please use data preparation function data_prep() first")
  }
  
  ## calculate appropriate transformations
  if (paired) {
    res <- tsf_paired(otu.data$otu.props, otu.data$otu.clr, otu.data$metadata, norm = norm)
  } else {
    res <- tsf_long(otu.data$otu.props, otu.data$otu.clr, otu.data$metadata, norm = norm)
    if (length(unique(table(otu.data$metadata$time))) != 1) { balanced = FALSE } else { balanced = TRUE }
  }
  if (paired) { type = "paired"
  } else if (balanced) { type = "balanced longitudinal"
  } else {
    type = "unbalanced longitudinal"
    warning("WARNING: this transformation is not recommended for strongly unbalanced designs!")
  }
  
  ## return
  return(list(dat.binary = res$dat.binary, dat.quant.prop = res$dat.quant.prop,
              dat.quant.clr = res$dat.quant.clr, avg.prop = res$avg.prop, type = type))
}


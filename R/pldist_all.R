#' pldist.all
#'
#' Function that calculates multiple paired and longitudinal ecological distance/dissimilarity 
#' matrices. Includes qualitative and quantitative versions of Bray-Curtis, Jaccard, Kulczynski, 
#' Gower, and unweighted and generalized UniFrac distances/dissimilarities. UniFrac-based 
#' metrics are based in part on GUniFrac (Jun Chen & Hongzhe Li (2012)). Both quantitative and 
#' qualitative versions of each requested metric are returned 
#'
#' @param otus OTU count or frequency table, containing one row per sample and one column per OTU. 
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier (if using longitudinal distances, this must be numeric or 
#'     convertable to numeric).
#' @param paired Logical indicating whether to use the paired version of the metric (TRUE) or the 
#'     longitudinal version (FALSE). Paired analyis is only possible when there are exactly 2 
#'     unique time points/identifiers for each subject or pair. 
#' @param method Desired distance metric(s). Include a vector with any combination of braycurtis,
#'     jaccard, kulczynski, gower, and unifrac, or any unambiguous abbreviation thereof. 
#' @param tree Rooted phylogenetic tree of R class "phylo". Default NULL; only needed for 
#'     UniFrac family distances. 
#' @param gam Parameter controlling weight on abundant lineages for UniFrac family distances. The 
#'     same weight is used within a subject as between subjects. Default (0, 0.5, 1); only needed for 
#'     UniFrac family distances. 
#' @return Returns a list containing all n x n distance (or dissimilarity) matrices requested, 
#' with both quantitative and qualitative versions of the metric, named as "D_metric_quant" or 
#' "D_metric_qual".
#'     
#'     
pldist_all <- function(otus, metadata, paired = FALSE, method = c("b", "g", "j", "k", "u"), 
                   tree = NULL, gam = c(0, 0.5, 1)) {
  ## Find desired method 
  method.opts = c("braycurtis", "jaccard", "kulczynski", "gower", "unifrac")
  this.method = pmatch(trimws(tolower(method)), method.opts, nomatch = NA)
  if (all(is.na(this.method))) {
    stop("Method does not match any expected methods. Please see list of options in documentation.")
  } 
  method = sort(method.opts[this.method])  # unifrac is last (if present) 
  
  okdat <- check_input(otus = otus, metadata = metadata, paired = paired)
  otus <- okdat$otus 
  metadata <- okdat$metadata 
  remove(okdat) 
  
  ## Calculate distances/dissimilarities 
  Ds <- list() 
  
  for (mm in 1:length(method)) {
    if (method[mm] != "unifrac") {
      ## Calculate transformed data and apply distance (all except UniFrac) 
      tsf.res <- pltransform(otus = otus, metadata = metadata, paired = paired, check.input = FALSE)
      this.D.bin <- switch(method[mm], 
                           braycurtis = braycurtis(tsf.res, binary = TRUE), 
                           jaccard = jaccard(tsf.res, paired = paired, binary = TRUE), 
                           kulczynski = kulczynski(tsf.res, paired = paired, binary = TRUE), 
                           gower = gower(tsf.res, binary = TRUE))
      this.D.quant <- switch(method[mm], 
                             braycurtis = braycurtis(tsf.res, binary = FALSE), 
                             jaccard = jaccard(tsf.res, paired = paired, binary = FALSE), 
                             kulczynski = kulczynski(tsf.res, paired = paired, binary = FALSE), 
                             gower = gower(tsf.res, binary = FALSE))      
      Ds[[(2*mm - 1)]] <- this.D.bin 
      Ds[[2*mm]] <- this.D.quant 
      names(Ds)[(2*mm-1):(2*mm)] <- paste("D", method[mm], c("qual", "quant"), sep = "_")
    } else {
      ## Calculate paired/longitudinal UniFrac dissimilarities 
      if (is.null(tree)) stop("Tree is required for UniFrac family metrics.")
      if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!") 
      this.D <- LUniFrac(otu.tab = otus, metadata = metadata, tree = tree, gam = gam, paired = paired, check.input = FALSE)
      
      Ds[[(2*mm - 1)]] <- this.D[,,"d_UW"]
      for (i in 1:length(gam)) {
        Ds[[(2*mm + i - 1)]] <- this.D[,, paste("d", gam[i], sep = "_")] 
      }
      names(Ds)[(2*mm-1)] <- "D_UW"
      names(Ds)[(2*mm):(2*mm + length(gam) - 1)] <- paste("D", gam, sep = "_")
    } 
  }
  
  return(Ds)
}

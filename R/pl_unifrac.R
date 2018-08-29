#' LUniFrac 
#'
#' Longitudinal UniFrac distances for comparing changes in
#' microbial communities across 2 time points.
#'
#' Based in part on Jun Chen & Hongzhe Li (2012), GUniFrac.
#'
#' Computes difference between time points and then calculates
#' difference of these differences, resulting in a dissimilarity
#' matrix that can be used in a variety of downstream 
#' distance-based analyses.
#'
#' @param otu.tab OTU count table, containing 2*n rows (samples) and q columns (OTUs)
#' @param tree Rooted phylogenetic tree of R class "phylo"
#' @param gam Parameter controlling weight on abundant lineages. The same weight is used within a subjects as between subjects.
#' @param metadata Data frame with three columns: subject identifiers (n unique values), 
#'     sample identifiers (must match row names of otu.tab), 
#'     and time or group indicator (numeric variable, or factor with levels such that as.numeric returns 
#'     the desired ordering). Column names should be subjID, sampID, time. 
#' @param paired Logical indicating whether to use the paired (TRUE) or longitudinal (FALSE) transformation. 
#' @return Returns a (K+1) dimensional array containing the longitudinal UniFrac dissimilarities 
#'    with the K specified gamma values plus the unweighted distance. The unweighted dissimilarity 
#'    matrix may be accessed by result[,,"d_UW"], and the generalized dissimilarities by result[,,"d_G"] 
#'    where G is the particular choice of gamma.
#' 
#' @export
#' 
LUniFrac <- function(otu.tab, tree, gam = c(0, 0.5, 1), metadata, paired) {
  n <- nrow(otu.tab)
  
  # Check OTU name consistency
  if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU names
         in the OTU table and the tree should match." )
  }
  
  # Get the subtree if tree contains more OTUs
  absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
  if (length(absent) != 0) {
    tree <- drop.tip(tree, absent)
    warning("The tree has more OTU than the OTU table!")
  }
  
  # Reorder the otu.tab matrix if the OTU orders are different
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]
  
  ntip <- length(tip.label)
  nbr <- nrow(tree$edge)        # number of branches = 2*(ntip - 1)
  edge <- tree$edge             # edges entering a node (1 through (ntip - 1))
  edge2 <- edge[, 2]            # edges leaving a node (1 through 2*(ntip -1))
  br.len <- tree$edge.length    # branch lengths, corresponds to edge2
  
  #  Accumulate OTU proportions up the tree
  ## Note: columns are samples, rows are branches 
  cum <- matrix(0, nbr, n)							# Branch abundance matrix (n = nsamp = nsubj * ntimes in this usage) 
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
    node <- edge[tip.loc, 1]						# Assume the direction of edge
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }
  colnames(cum) = rownames(otu.tab)
  
  ### Step 1: calculate within-subject distance data
  if (paired) {
    tsf.dat <- pltransform(otus = t(cum), metadata = metadata, paired = TRUE, check.input = FALSE)$tsf.data 
  } else {
    tsf.dat <- pltransform(otus = t(cum), metadata = metadata, paired = FALSE, check.input = FALSE)$tsf.data 
  }
  cum.avg <- t(tsf.dat$avg.prop)
  cum.unw <- t(tsf.dat$dat.binary)
  cum.gen <- t(tsf.dat$dat.quant) 
  
  # Construct the returning array
  # d_UW: unweighted
  dimname3 <- c(paste("d", gam, sep="_"), "d_UW")
  lunifracs <- array(NA, c(ncol(cum.avg), ncol(cum.avg), length(gam) + 1),
                     dimnames=list(colnames(cum.avg), colnames(cum.avg), dimname3))
  for (i in 1:(length(gam)+1)){
    for (j in 1:ncol(cum.avg)){
      lunifracs[j, j, i] <- 0
    }
  }
  
  ### Step 2: calculate distances based on within-subject summaries
  for (i in 2:ncol(cum.avg)) {
    for (j in 1:(i-1)) {
      d1 <- cum.gen[, i]
      d2 <- cum.gen[, j]
      avg1 <- cum.avg[, i]
      avg2 <- cum.avg[, j]
      
      ind <- which((abs(d1) + abs(d2)) != 0)
      d1 <- d1[ind]
      d2 <- d2[ind]
      diff <- abs(d2 - d1)
      br.len2 <- br.len[ind]
      
      # Generalized LUniFrac dissimilarity
      for(k in 1:length(gam)){
        w <- br.len * (avg1 + avg2)^gam[k]
        lunifracs[i, j, k] = lunifracs[j, i, k] = sum(diff * w[ind]) / sum(w)
      }
      
      #	Unweighted LUniFrac Distance
      d1 <- cum.unw[, i]
      d2 <- cum.unw[, j]
      diff <- abs(d2 - d1) 
      
      # only branches with some change contribute
      ind <- which((abs(cum.unw[, i]) + abs(cum.unw[,j])) != 0)
      if (length(ind) > 0) {
        diff <- diff[ind]
        br.len2 <- br.len[ind]
        lunifracs[i, j, (k + 1)] = lunifracs[j, i, (k + 1)] = sum(br.len2*diff) / sum(br.len)
      } else {
        lunifracs[i, j, (k + 1)] = lunifracs[j, i, (k + 1)] = 0
      }
    }
  }
  return(lunifracs)
}
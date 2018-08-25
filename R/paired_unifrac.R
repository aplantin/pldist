#' PUniFrac
#'
#' Paired UniFrac distances for comparing changes in
#' microbial communities across 2 groups or time points.
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
#' @param gam Parameter controlling weight on abundant lineages. The same weight is used within a subject as between subjects.
#' @param metadata Data frame with three columns: subject identifiers (n unique values, column name "subjID"), 
#'     sample identifiers (must match row names of otu.tab, column name "sampID"), 
#'     and time point or group identifier (variable with two unique levels, column name "time").
#' @return Returns a (K+1) dimensional array containing the longitudinal UniFrac dissimilarities 
#'     with the K specified gamma values plus the unweighted distance. The unweighted dissimilarity 
#'     matrix may be accessed by result[,,"d_UW"], and the generalized dissimilarities by result[,,"d_G"] 
#'     where G is the particular choice of gamma.
#' @importFrom ape is.rooted drop.tip rtree
#' 
#' @export
#'
PUniFrac <- function(otu.tab, tree, gam = c(0, 0.5, 1), metadata) {
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
  cum <- matrix(0, nbr, n)							# Branch abundance matrix
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
  
  ### Step 1: calculate within-subject distance data
  metadata[,3] <- as.numeric(as.factor(metadata[,3]))
  cum.t1 <- cum[, which(metadata[,3] == 1)]
  colnames(cum.t1) <- metadata[,1][which(metadata[,3] == 1)]
  
  cum.t2 <- cum[, which(metadata[,3] == 2)]
  colnames(cum.t2) <- metadata[,1][which(metadata[,3] == 2)]
  cum.t2 <- cum.t2[, colnames(cum.t1)]
  
  if (any(colnames(cum.t1) != colnames(cum.t2))) {
    stop("Same set of samples is not present at both time points!")
  }
  
  cum.avg <- Reduce("+", list(cum.t1, cum.t2))/2
  
  # Construct the returning array
  # d_UW: unweighted
  dimname3 <- c(paste("d", gam, sep="_"), "d_UW")
  lunifracs <- array(NA, c(ncol(cum.t1), ncol(cum.t1), length(gam) + 1),
                     dimnames=list(colnames(cum.t1), colnames(cum.t1), dimname3))
  for (i in 1:(length(gam)+1)){
    for (j in 1:ncol(cum.t2)){
      lunifracs[j, j, i] <- 0
    }
  }
  
  # Calculate within-subject distances (for generalized/weighted)
  cum.diff <- matrix(0, nrow = nrow(cum), ncol = ncol(cum.t2))
  for (i in 1:ntip) {
    cum.diff[i,] <- (cum.t2[i,] - cum.t1[i,]) / (cum.t2[i,] + cum.t1[i,])
  }
  cum.diff[is.na(cum.diff)] <- 0
  
  # Calculate within-subject distances (for unweighted)
  cum.diff.uw <- matrix(0, nrow = nrow(cum), ncol = ncol(cum.t2))
  for (i in 1:ntip) {
    ## this is the change in V2
    cum.diff.uw[i,] <- as.numeric(cum.t2[i,] > 0) - as.numeric(cum.t1[i,] > 0)
  }
  
  ### Step 2: calculate distances based on within-subject summaries
  for (i in 2:ncol(cum.diff)) {
    for (j in 1:(i-1)) {
      d1 <- cum.diff[, i]
      d2 <- cum.diff[, j]
      avg1 <- cum.avg[, i]
      avg2 <- cum.avg[, j]
      
      ind <- which((abs(d1) + abs(d2)) != 0)
      d1 <- d1[ind]
      d2 <- d2[ind]
      avg1 <- avg1[ind]
      avg2 <- avg2[ind]
      br.len2 <- br.len[ind]
      
      diff <- abs(d2 - d1)/2 
      
      # Generalized LUniFrac dissimilarity
      for(k in 1:length(gam)){
        w <- br.len2 * (avg1 + avg2)^gam[k]
        lunifracs[i, j, k] = lunifracs[j, i, k] = sum(diff * w) / sum(w)
      }
      
      #	Unweighted LUniFrac Distance
      d1 <- cum.diff.uw[, i]
      d2 <- cum.diff.uw[, j]
      diff <- abs(d2 - d1)/2
      
      # only branches with some change contribute
      ind <- which((abs(cum.diff.uw[, i]) + abs(cum.diff.uw[,j])) != 0)
      if (length(ind) > 0) {
        diff <- diff[ind]
        br.len2 <- br.len[ind]
        lunifracs[i, j, (k + 1)] = lunifracs[j, i, (k + 1)] = sum(br.len2*diff) / sum(br.len2)
      } else {
        lunifracs[i, j, (k + 1)] = lunifracs[j, i, (k + 1)] = 0
      }
    }
  }
  return(lunifracs)
}
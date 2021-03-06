#' clr_LUniFrac 
#'
#' Longitudinal UniFrac distances for comparing changes in
#' microbial communities across 2 time points, using CLR-transformed data.
#'
#' Based in large part on Jun Chen & Hongzhe Li (2012), GUniFrac. 
#' With reference to the GitHub account of user ruthgrace, repository 
#' CLRUniFrac (this method is not associated with a publication). 
#'
#' Computes difference between time points and then calculates
#' difference of these differences, resulting in a dissimilarity
#' matrix that can be used in a variety of downstream 
#' distance-based analyses.
#'
#' @param otu.tab OTU count table, containing 2*n rows (samples) and q columns (OTUs)
#' @param metadata Data frame with three columns: subject identifiers (n unique values), 
#'    sample identifiers (must match row names of otu.tab), 
#'    and time or group indicator (numeric variable, or factor with levels such that as.numeric returns 
#'    the desired ordering). Column names should be subjID, sampID, time. 
#' @param tree Rooted phylogenetic tree of R class "phylo"
#' @param gam Parameter controlling weighting factor for average taxon abundance. 
#' @param paired Logical indicating whether to use the paired (TRUE) or longitudinal (FALSE) transformation. 
#' @param pseudocount Pseudocount to be added to all values in OTU matrix prior to CLR transformation. 
#'     Default NULL. If NULL, then 0.5 is added if entries are counts, or min(1e-06, 0.5*min(nonzero p)) 
#'     if entries are proportions. If all entries are nonzero, nothing is added. 
#'     
#' @importFrom phytools getDescendants 
#' 
#' @return Returns a (K+1) dimensional array containing the longitudinal UniFrac dissimilarities 
#'    with the K specified gamma values plus the unweighted distance. The unweighted dissimilarity 
#'    matrix may be accessed by result[,,"d_UW"], and the generalized dissimilarities by result[,,"d_G"] 
#'    where G is the particular choice of gamma.
#'    
#' @examples
#' data("bal.long.otus")
#' data("bal.long.meta")
#' data("sim.tree")
#' D2.unifrac <- clr_LUniFrac(otu.tab = bal.long.otus, metadata = bal.long.meta, 
#'     tree = sim.tree, gam = c(0, 0.5, 1), paired = FALSE)
#' D2.unifrac[, , "d_1"]   # gamma = 1 (quantitative longitudinal transformation)
#' D2.unifrac[, , "d_UW"]  # unweighted LUniFrac (qualitative/binary longitudinal transf.)
#' 
#' @export
#' 
clr_LUniFrac <- function(otu.tab, metadata, tree, gam = c(0, 0.5, 1), paired, pseudocount = NULL) {
  # check data 
  temp <- data_prep(otu.tab, metadata, paired)  # just for data checking, not transformations
  
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
  
  # Add pseudocount and re-close to proportions 
  otu.psct = otu.tab
  if (any(otu.psct == 0)) {
    if (is.null(pseudocount)) {
      subjsums <- apply(otu.psct, 1, FUN = function(x) sum(x))
      if (all(subjsums == 1)) {  # data are proportions 
        otu.psct <- otu.psct + min(1e-06, otu.psct[otu.psct != 0])
      } else { otu.psct <- otu.psct + 0.5 }
    } else { otu.psct <- otu.psct + pseudocount }
  }
  otu.psct <- counts2props(otu.psct)
  otu.tab <- counts2props(otu.tab) 
 
  # Calculate tip-level CLR transformed data and geometric mean for each sample 
  samp.gm <- apply(otu.psct, 1, FUN = function(x) mean(log(x))) 
  otu.tab.clr <- t(apply(otu.psct, 1, FUN = function(x) log(x) - mean(log(x))))
  nsamp <- nrow(otu.psct) 

  # Summarize tree 
  ntip <- length(tip.label)     # number of OTUs (# tips on tree = # columns)
  nbr <- nrow(tree$edge)        # number of branches = 2*(ntip - 1)
  edge <- tree$edge             # edges entering a node (1 through (ntip - 1))
  edge2 <- edge[, 2]            # edges leaving a node (1 through 2*ntip -1)
  nn <- 2*ntip - 1 
  br.len <- tree$edge.length    # branch lengths, corresponds to (leaving) edge2 node 
  
  # Cumulative proportions and CLR up the tree 
  ## Note: columns are samples, rows are branches (transpose of normal OTU matrix) 
  cum.psct <- matrix(0, nbr, nsamp)	# Branch abundance matrix (nsamp = nsubj * ntimes)
  cum.orig <- matrix(0, nbr, nsamp) 
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)    # the row of `edge` corresponding to branch leaving this tip 
    cum.psct[tip.loc, ] <- cum.psct[tip.loc, ] + otu.psct[, i]
    cum.orig[tip.loc, ] <- cum.orig[tip.loc, ] + otu.tab[, i] 
    node <- edge[tip.loc, 1]						# Assume the direction of edge
    node.loc <- which(edge2 == node)    # next-level-up node (the other end of this edge/branch)
    while (length(node.loc)) {
      cum.psct[node.loc, ] <- cum.psct[node.loc, ] + otu.psct[, i]
      cum.orig[node.loc, ] <- cum.orig[node.loc, ] + otu.tab[, i]
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }
  colnames(cum.psct) = rownames(otu.psct)
  colnames(cum.orig) = rownames(otu.tab) 
  rownames(cum.psct) = rownames(cum.orig) = edge2
  
  # cumulative CLR-transformed 
  clr.cum = matrix(0, nbr, nsamp)
  colnames(clr.cum) = rownames(otu.psct)
  rownames(clr.cum) = edge2
  for (i in 1:nbr) {
    this.children <- getDescendants(tree, rownames(cum.psct)[i])
    if (length(this.children) > 0) {
      this.idx <- unique(c(rownames(cum.psct)[i], (1:ntip)[-which((1:ntip) %in% this.children)]))
    } else {
      this.idx <- unique(c(rownames(cum.psct)[i], (1:ntip)))
    }
    if (length(this.idx) > 1) {
      this.gm <- apply(cum.psct[this.idx, ], 2, FUN = function(x) mean(log(x)))  # GM if all under given node are this OTU 
      clr.cum[i, ] <- log(cum.psct[i, ]) - this.gm
    } else {  # should only be root node 
      clr.cum[i, ] <- 0  
    }
  }

  # given cumulative CLR, calculate within-subject distances 
  if (paired) {
    tsf.dat <- pltransform(otu.data = list(otu.props = t(cum.orig), otu.clr = t(clr.cum), metadata = metadata), paired = TRUE) 
  } else {
    tsf.dat <- pltransform(otu.data = list(otu.props = t(cum.orig), otu.clr = t(clr.cum), metadata = metadata), paired = FALSE) 
  }
  cum.avg <- t(tsf.dat$avg.prop)
  cum.unw <- t(tsf.dat$dat.binary)
  cum.gen <- t(tsf.dat$dat.quant.clr) 
  
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


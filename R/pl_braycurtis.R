#' Paired or longitudinal Bray-Curtis distances 
#' 
#' The distances are calculated as follows, where d_k^X is the within-subject 
#'     measure of change appropriate to the setting (paired/longitudinal and 
#'     quantitative/qualitative), as described in the full package documentation 
#'     and vignette. 
#' \eqn{D_{AB} = (1/m) * \sum_k |d_k^A - d_k^B|}
#'
#' @param tsf.data Transformed OTU table and metadata (from function pl.transform)
#' @param binary Logical indicating whether to use the binary version of the distance 
#' @return Returns an n x n distance matrix. 
#'
#' @export
#' 
braycurtis <- function(tsf.data, binary) {
  if (binary) { dat = tsf.data$dat.binary
  } else { dat = tsf.data$dat.quant }
  
  n = nrow(dat); m = ncol(dat) 
  out.D <- matrix(0, n, n) 
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      out.D[i, j] = out.D[j, i] = 1/m * sum(abs(dat[i, ] - dat[j, ]))
    }
  }
  
  return(out.D) 
}
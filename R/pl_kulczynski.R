#' Paired or longitudinal Kulczynski distances 
#'
#' The distances are calculated as follows, where d_k^X is the within-subject 
#'     measure of change appropriate to the setting (paired/longitudinal and 
#'     quantitative/qualitative), as described in the full package documentation 
#'     and vignette. 
#'     
#' Paired, qualitative: \eqn{D_{AB} = 1 - (1/m) \sum_k I[d_k^A = d_k^B] I[d_k^A \neq 0]} 
#' Paired, quantitative: \eqn{D_{AB} = 1 - (2/m) \sum_k \min(|d_k^A|, |d_k^B|) I[sgn(d_k^A) = sgn(d_k^B)]}
#' Longitudinal: \eqn{D_{AB} = 1 - (1/m) * \sum_k \min(d_k^A, d_k^B)}
#'    
#' @param tsf.data Transformed OTU table and metadata (from function pl.transform)
#' @param paired Logical indicating whether paired analysis is desired 
#' @param binary Logical indicating whether to use the binary version of the distance 
#' @return Returns an n x n distance matrix. 
#'
#' @export
#' 
kulczynski <- function(tsf.data, paired, binary) {
  if (binary) { dat = tsf.data$dat.binary
  } else { dat = tsf.data$dat.quant }
  
  n = nrow(dat); m = ncol(dat) 
  out.D <- matrix(0, n, n) 
  
  if (paired) {
    if (binary) {
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          out.D[i, j] = out.D[j, i] = 1 - sum(dat[i, ] == dat[j, ] & dat[i, ] != 0)/m 
        }
      }
    } else {   # paired but not binary 
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          num = sum(pmin(abs(dat[i,]), abs(dat[j,])) * flexsign(dat[i, ], dat[j, ]))
          out.D[i, j] = out.D[j, i] = 1 - 2*num/m
        }
      }
    }
  } else {  # not paired (longitudinal) 
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        out.D[i,j] = out.D[j,i] = 1 - sum(pmin(dat[i,],dat[j,])) / m
      }
    }
  }
  return(out.D) 
}
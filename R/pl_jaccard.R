#' flexsign
#' 
#' Sign function that considers 0 both positive and negative. Returns 1 if the two numbers are the same sign, 0 otherwise. Vectorized (compares vectors elementwise). 
#'
#' @param v1 First vector
#' @param v2 Second vector
#' @return Returns an n x n distance matrix. 
#'
#' @export
#' 
flexsign <- function(v1, v2) {
  return(as.numeric( (v1 >= 0 & v2 >= 0) | (v1 <= 0 & v2 <= 0) ))
}


#' Paired or longitudinal Jaccard distances 
#' 
#' The distances are calculated as follows, where d_k^X is the within-subject 
#'     measure of change appropriate to the setting (paired/longitudinal and 
#'     quantitative/qualitative), as described in the full package documentation 
#'     and vignette. 
#' Paired, qualitative: \eqn{D_{AB} = 1 - {\sum_k I(d_k^A = d_k^B) I(d_k^A \neq 0)}/{\sum_k [I(d_k^A \neq 0) + I(d_k^B \neq 0)]}}
#' Paired, quantitative: \eqn{D_{AB} = 1 - {\sum_k \min(|d_k^A|, |d_k^B|) \, I(sgn(d_k^A) = sgn(d_k^B))}/{\sum_k \max(|d_k^A|, |d_k^b|)}}
#' Longitudinal: \eqn{D_{AB} = 1 - (\sum_k \min(d_k^A, d_k^B))/(\sum_k \max(d_k^A, d_k^B))}
#'
#' @param tsf.data Transformed OTU table and metadata (from function pl.transform)
#' @param paired Logical indicating whether paired analysis is desired 
#' @param binary Logical indicating whether to use the binary version of the distance 
#' @return Returns an n x n distance matrix. 
#'
#' @export
#' 
jaccard <- function(tsf.data, paired, binary) {
  if (binary) { dat = tsf.data$dat.binary
  } else { dat = tsf.data$dat.quant }
  
  n = nrow(dat); m = ncol(dat) 
  out.D <- matrix(0, n, n) 
  
  if (paired) {
    if (binary) {
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          num = sum(as.numeric(dat[i, ] == dat[j, ] & dat[i, ] != 0))
          denom = sum(as.numeric(dat[i, ] != 0)) + sum(as.numeric(dat[j, ] != 0))
          if (denom != 0) {
            out.D[i,j] = out.D[j,i] = 1 - num/denom
          } else {
            out.D[i,j] = out.D[j,i] = 0
          }
          
        }
      }
    } else { # paired, not binary 
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          idx = which(pmax(abs(dat[i,]), abs(dat[j,])) != 0)
          num = sum(pmin(abs(dat[i,idx]), abs(dat[j,idx])) * flexsign(dat[i,idx], dat[j,idx]))
          denom = sum(pmax(abs(dat[i,idx]), abs(dat[j,idx])))
          out.D[i, j] = out.D[j, i] = 1 - num/denom
        }
      }
    }
  } else {  # not paired 
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        idx = which(pmax(dat[i,], dat[j,]) != 0)
        out.D[i,j] = out.D[j,i] = 1 - sum(pmin(dat[i,idx],dat[j,idx])) / sum(pmax(dat[i,idx],dat[j,idx]))
      }
    }
  }
  return(out.D) 
}
#' Function to recreate the upper triangular matrix
#'
#' @param D the tropical interpretation of a tree as a vector
#' @param n a numeric value indicating the the number of leaves of the tree
#' @param tips a vector of leaf names
#' @return A cophenetic distance matrix, where each entry is the distance
#' between any two leaves of the tree
#' @keywords internal
make.matrix <- function(D, n, tips){
  dd <- matrix(rep(0, n*n), nrow=n, ncol=n, byrow=TRUE)
  count <- 1
  for(i in 1:(n-1))
    for(j in (i + 1):n){
      dd[i, j] <- D[count]
      dd[j, i] <- dd[i, j]
      count <- count+1
    }
  mymatrix <- matrix(dd, nrow=n, ncol=n, byrow=TRUE, dimnames=list(tips,tips))
  return(mymatrix)
}

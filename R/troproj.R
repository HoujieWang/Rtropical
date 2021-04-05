#' Projection on Tropical Polytope
#'
#' \code{troproj} computes a vector projection on a given tropical polytope.
#'
#' @param x a data vector.
#' @param tconv a data matrix, of dimension e x s; e is the dimension of the tropical space;
#' s is the dimension of the polytope. each column is a basis of tropical polytope.
#'
#' @return A projected vector on the given tropical polytope.
#'
#' @author Houjie Wang
#'
#' @references Page, R., Yoshida, R. & Zhang L.
#' \emph{Tropical principal component analysis on the space of phylogenetic trees.
#' J. Bioinform., Volume 36, Issue 17, 4590–4598 (2020).}
#' \url{https://doi.org/10.1093/bioinformatics/btaa564}
#'
#' @examples
#'
#' # Generate a tropical polytope consisting of three trees each with 5 leaves
#' library(ape)
#' pltp <- sapply(1: 3, function(i){vec.fun(rcoal(5))})
#' # Generate an observation and vectorize it
#' tree <- rcoal(5)
#' tree_vec <- vec.fun(tree)
#' troproj(tree_vec, pltp)
#'
#' @export
#' @export troproj
#'
troproj <- function(x, tconv){
  if(is.null(dim(tconv))){
    lambda <- min(x - tconv)
    pi_D <- c(t(lambda + t(tconv)))
  }else{
    lambda <- apply(x - tconv, 2, min)#D_s by row
    pi_D <- apply(t(lambda + t(tconv)),1,max)
  }
  return(pi_D)
}

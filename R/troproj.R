#' Projection on Tropical Polytope
#'
#' \code{troproj} computes a projection onto a given tropical polytope.
#'
#' @param x x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param tconv a data matrix, of dimension nvars x s; each column is a basis of tropical polytope;
#' s is the dimension of the polytope.
#'
#' @return The projected vector on the given tropical polytope.
#'
#' @author Qiwen Kang and Houjie Wang
#'
#' Maintainer: Houjie Wang \email{wanghoujie6688@@gmail.com}
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
troproj <- function(x,tconv){
  if(is.null(dim(x))){
    lambda <- min(tconv - x)
    pi_D <- c(t(lambda + t(x)))
  }else{
    lambda <- apply(tconv - x, 2, min)
    pi_D <- apply(t(lambda + t(x)),1,max)
  }
  return(pi_D)
}

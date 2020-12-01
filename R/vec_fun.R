#' Vectorize a tree
#'
#' @importFrom ape cophenetic.phylo
#' @param x A object of class "phylo" in R
#' @export
#' @examples
#' vec_fun(tree)
#'
vec_fun<-function(x){
  x = cophenetic.phylo(x)
  x[lower.tri(x)]
}

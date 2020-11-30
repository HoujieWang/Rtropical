#' Vectorize a tree
#'
#' @importFrom ape cophenetic.phylo
#' @param x A tree object in R as list
#' @export
#' @examples
#' vec_fun(tree)
#'
vec_fun<-function(x){
  x = cophenetic.phylo(x)
  x[lower.tri(x)]
}

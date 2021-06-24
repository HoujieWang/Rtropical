#' Vectorize a tree
#'
#' @importFrom ape cophenetic.phylo
#' @param phy A object of class "phylo" in R
#' @export
#'
vec.fun <- function(phy) {
  phy <- cophenetic.phylo(phy)
  phy[lower.tri(phy)]
}

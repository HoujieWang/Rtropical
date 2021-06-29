#' Vectorize a tree
#'
#' @importFrom ape cophenetic.phylo
#' @param phy A object of class phylo
#' @export
#'
vec.fun <- function(phy) {
  cophe <- cophenetic.phylo(phy)
  labs <- outer(rownames(cophe), colnames(cophe), paste, sep="--")
  phyvec <- cophe[lower.tri(cophe)]
  names(phyvec) <- labs[lower.tri(cophe)]
  phyvec
}

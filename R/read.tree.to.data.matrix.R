#' Read trees in two categories into a data matrix
#'
#' @importFrom ape read.tree
#' @param data.file1 A set of trees in the format of Newick strings from one category.
#' @param data.file2 a set of trees in the format of Newick strings from another category.
#' @export
#'
read.tree.to.data.matrix <- function(data.file1, data.file2){
  G1 <- read.tree(data.file1)
  G2 <- read.tree(data.file2)

  n <- length(G1[[1]]$tip.label)
  to <- G1[[1]]$tip.label
  N1 <- length(G1)
  N2 <- length(G2)

  distVec_all1 <- multiPhylo.to.data.matrix(G1, to)
  distVec_all2 <- multiPhylo.to.data.matrix(G2, to)
  rownames(distVec_all1) <- NULL
  rownames(distVec_all2) <- NULL

  class <- as.factor(c(rep(1, N1), rep(2, N2)))
  D <- cbind(class, rbind(distVec_all1, distVec_all2))
  return(D)
}

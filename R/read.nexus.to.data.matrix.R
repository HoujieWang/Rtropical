#' Read trees in two categories into a data matrix
#'
#' @param data.file1 A data set with trees from one category.
#' @param data.file2 A data set with trees from the other category.
#' @export
#' @examples
#' read.nexus.to.data.matrix(data.file1, data.file2)
#'
read.nexus.to.data.matrix <- function(data.file1, data.file2){
  # data.file1 = "genetree_S1_r025.dat"
  # data.file2 = "genetree_S2_r025.dat"
  G1 <- read.nexus(data.file1)
  G2 <- read.nexus(data.file2)

  n <- length(G1[[1]]$tip.label)
  to <- G1[[1]]$tip.label
  N1 <- length(G1)
  N2 <- length(G2)

  distVec_all1 <- distMat(G1, tipOrder = to)
  distVec_all2 <- distMat(G2, tipOrder = to)

  D_all1 <- matrix(unlist(distVec_all1), ncol=N1)
  D_all2 <- matrix(unlist(distVec_all2), ncol=N2)
  class1 <- rep(1, dim(D_all1)[2])
  class2 <- rep(2, dim(D_all2)[2])
  D_all1 <- rbind(class1, D_all1)
  D_all2 <- rbind(class2, D_all2)

  D <- rbind(t(D_all1), t(D_all2))

  return(D)
}

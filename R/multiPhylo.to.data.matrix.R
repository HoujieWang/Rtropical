#' Convert an object in class "multiPhylo" into a data matrix
#'
#' @importFrom parallel parLapply
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel setDefaultCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom ape root
#' @importFrom ape chronos
#' @param phy An object in class "multiPhylo"
#' @param tipOrder Order of leaf names
#' @return A data matrix with the first column of the categories
#' @export
#'
multiPhylo.to.data.matrix <- function(phy, tipOrder){
  cl <- makeCluster(2)
  setDefaultCluster(cl)
  clusterExport(cl, c("chronos", "vec.fun"), envir = environment())
  trees_root <- root(phy, outgroup = tipOrder[1],resolve.root=TRUE)
  chronotrees <- parLapply(cl, trees_root, chronos)
  distVec_all <- parLapply(cl, chronotrees, vec.fun)
  stopCluster(cl)
  return(do.call("rbind", distVec_all))
}

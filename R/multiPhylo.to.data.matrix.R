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
#' @param trees An object in class "multiPhylo" containing trees
#' @param tipOrder Order of leaf names
#' @export
#' @examples
#' multiPhylo.to.data.matrix(trees, tipOrder)
#'
multiPhylo.to.data.matrix <- function(trees, tipOrder){
  cl <- makeCluster(2)
  setDefaultCluster(cl)
  clusterExport(cl, c("chronos", "vec.fun"))
  if(class(trees) == "phylo"){
    trees = list(trees)
    class(trees) = "multiPhylo"
  }
  tipOrder = trees[[1]]$tip.label
  trees_root <- root(trees, outgroup = tipOrder[1],resolve.root=TRUE)
  chronotrees <- parLapply(cl, trees_root, chronos)
  distVec_all <- parLapply(cl, chronotrees, vec.fun)
  stopCluster(cl)
  return(do.call("rbind", distVec_all))
}

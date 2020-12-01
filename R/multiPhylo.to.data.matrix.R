#' Convert an object in class "multiPhylo" into a data matrix
#'
#' @importFrom parallel mclapply
#' @param trees An object in class "multiPhylo" containing trees
#' @export
#' @examples
#' multiPhylo.to.data.matrix(trees)
#'
multiPhylo.to.data.matrix <- function(trees){
  if(class(trees) == "phylo"){
    trees = list(trees)
    class(trees) = "multiPhylo"
  }
  tipOrder = trees[[1]]$tip.label
  trees_root <- root(trees, outgroup = tipOrder[1],resolve.root=TRUE)
  chronotrees <- mclapply(trees_root, chronos)
  distVec_all <- mclapply(chronotrees, vec.fun)
  return(do.call("rbind", distVec_all))
}

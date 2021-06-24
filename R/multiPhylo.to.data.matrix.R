#' Convert an object in class "multiPhylo" into a data matrix
#'
#' @importFrom parallel parLapply
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel setDefaultCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom ape root
#' @importFrom ape chronos
#' @param phy an object in class "multiPhylo"
#' @param tipOrder a numeric vector of order of leaf names
#' @param ncores a numeric value indicating the number of threads
#' utilized for multi-cored CPUs. (default: 2)
#' @return a data matrix with the first column of the categories
#' @export
#'
multiPhylo.to.data.matrix <- function(phy, tipOrder, ncores = 2) {
  if (Sys.info()["sysname"] == "Windows") {
    cl <- makeCluster(2)
    setDefaultCluster(cl)
    clusterExport(cl, c("chronos", "vec.fun"), envir = environment())
    trees_root <- root(phy, outgroup = tipOrder[1], resolve.root = TRUE)
    chronotrees <- parLapply(cl, trees_root, chronos)
    distVec_all <- parLapply(cl, chronotrees, vec.fun)
    stopCluster(cl)
  } else {
    trees_root <- root(phy, outgroup = tipOrder[1], resolve.root = TRUE)
    chronotrees <- mclapply(trees_root, chronos, mc.cores = ncores)
    distVec_all <- mclapply(chronotrees, vec.fun, mc.cores = ncores)
  }
  return(do.call("rbind", distVec_all))
}

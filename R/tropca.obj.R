#' Compute Tropical PCA Objective
#'
#' @importFrom parallel parLapply
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast rowMins
#'
#' @param pc a matrix of principle components
#' @param x_list a list of vectors
#' @param cl cluster for parallel computing
#' @export
#' @export tropca.obj
#'
tropca.obj <- function(pc, x_list, cl){
  proj <- parLapply(cl, x_list, troproj , tconv = pc)
  temp = do.call("rbind", x_list) - do.call("rbind", proj)
  sum(rowMaxs(temp, value = T) - rowMins(temp, value = T))
}

#' Projection on Tropical Linear Space
#'
#' troproj.linsp computes a vector projection on a given tropical linear space
#'
#' @importFrom Rfast rowSort
#'
#' @param x a data vector or a data matrix with each row an observation.
#' @param V a data matrix, of dimension s x e; e is the dimension of the tropical space;
#' s is the dimension of the linear space.
#'
#' @export
#' @export troproj.linsp
troproj.linsp = function(x, V){
  if(is.vector(x)) x = t(as.matrix(x))
  n = nrow(x); pcs = nrow(V); e = ncol(x);
  all_dets = array(NA, dim = c(e-2, e-1, e))
  all_combns = comboGeneral(1: e, pcs)
  all_combns_list = lapply(1: nrow(all_combns), function(i){all_combns[i, ]})
  all_dets[all_combns] = unlist(lapply(all_combns_list, function(i){tropdet(V[, i])}))

  data_proj = matrix(0, nrow = n, ncol = e)
  for (i in 1: e){
    # i = 1
    all_tau = comboGeneral(c(1: e)[-i], pcs - 1)
    all_tau = lapply(1: nrow(all_tau), function(i){all_tau[i, ]})
    temp2 = lapply(all_tau, function(tau){
      all_j = c(1: e)[-tau]
      temp_block = eachrow(x[, all_j], all_dets[rowSort(cbind(all_j, matrix(tau, ncol = 2, nrow = (e-pcs+1), byrow = T)))], "-")
      return(rowMins(temp_block, T) + tropdet(V[, c(tau, i)]))
    })
    temp2 = do.call(cbind, temp2)
    data_proj[, i] = rowMaxs(temp2, T)
  }
  data_proj
}

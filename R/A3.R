A3 <- function(data,tst_data,indice){
  library(lpSolve);library(gtools)
  count = 0
  k1 = indice[1]; k2 = indice[2] 
  n = nrow(data)/2; ntst = nrow(tst_data)
  p_tst = rep(0,ncol(tst_data)); q_tst = rep(0,ncol(tst_data))
  omega = rep(0, ncol(tst_data))
  obj = c(1, 0, 0, rep(-1, 4*n))
  conp = rbind(cbind(rep(1, n), matrix(rep(c(-1, 1), n), nrow = n, ncol = 2, byrow = T)), 
               cbind(rep(0, n), matrix(rep(c(-1, 1), n), nrow = n, ncol = 2, byrow = T)))
  conq = rbind(cbind(rep(1, n), matrix(rep(c(1, -1), n), nrow = n, ncol = 2, byrow = T)), 
               cbind(rep(0, n), matrix(rep(c(1, -1), n), nrow = n, ncol = 2, byrow = T)))
  con = cbind(rbind(conp, conq), diag(-1, nrow = 4*n, ncol = 4*n))
  rhs = c(rep(data[1: n, k1] - data[1: n, k2], 2), rep(data[-c(1: n), k2] - data[-c(1: n), k1], 2))
  dir=rep("<=",nrow(con))
  sol = lp("max", obj, con, dir, rhs)$solution
  omega[c(k1, k2)] = sol[2: 3]
  omega[-c(k1, k2)] = apply(data[, -c(k1, k2)], 2, function(x){min(c(min(data[1:n, k2] - x[1:n])+omega[k2], 
                                                                     min(data[(n+1):(2*n), k1] - x[(n+1):(2*n)])+omega[k1]))})
  tst_matrix = t(apply(tst_data, 1, function(x){x+omega}))
  tst_max = apply(tst_matrix,1,max)
  for (i in 1:ntst){
    if (tst_matrix[i,k2] == tst_max[i] | (tst_matrix[i,k1] != tst_max[i] & tst_matrix[i,k2] != tst_max[i])){
      q_tst = rbind(q_tst,tst_data[i,])
      if (i>(ntst/2)){
        count = count + 1
      }
    } else{
      p_tst = rbind(p_tst,tst_data[i,])
      if (i<=(ntst/2)){
        count = count + 1
      }
    }
  }
  return(list("omega " = omega,"best_accuracy" = count/ntst))
}
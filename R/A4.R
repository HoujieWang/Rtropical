A4 <- function(data,tst_data,indice){
  library(lpSolve);library(gtools)
  ip = indice[1]; j = indice[2]; iq = indice[3]
  n = nrow(data)/2; ntst = nrow(tst_data)/2
  obj = c(1, 0, 0, 0, rep(-1, 6*n))
  conp = rbind(cbind(rep(1, n), rep(-1, n), rep(0, n), rep(1, n)), cbind(rep(0, n), rep(-1, n), rep(0, n), rep(1, n)), cbind(rep(0, n), rep(0, n), rep(1, n), rep(-1, n)))
  conq = rbind(cbind(rep(1, n), rep(0, n), rep(-1, n), rep(1, n)), cbind(rep(0, n), rep(0, n), rep(-1, n), rep(1, n)), cbind(rep(0, n), rep(1, n), rep(0, n), rep(-1, n)))
  con = cbind(rbind(conp, conq), diag(-1, nrow = 6*n, ncol = 6*n))
  rhs = c(rep(data[1: n, ip] - data[1: n, j], 2), data[1: n, j] - data[1: n, iq], rep(data[-c(1: n), iq] - data[-c(1: n), j], 2), data[-c(1: n), j] - data[-c(1: n), ip])
  dir=rep("<=",nrow(con))
  omega = rep(0, ncol(tst_data))
  sol = lp("max", obj, con, dir, rhs)$solution
  omega[c(ip, iq, j)] = sol[2: 4]
  omega[-c(ip, iq, j)] = apply(data[, -c(ip, iq, j)], 2, function(x){min(data[, j] - x) + omega[j]})
  classification = apply(t(apply(tst_data, 1, function(x){x+omega})), 1, function(x){which(x == max(x))})
  PQ_com = matrix(c(1, 0, 0,
                    0, 1, 0,
                    0, 0, 1,
                    1, 1, 0,
                    1, 0, 1,
                    0, 1, 1,
                    1, 1, 1,
                    0, 0, 0), ncol = 3, byrow = T)
  ind_matrix = combinations(8, 4)
  accuracy = c()
  for (i in 1:nrow(ind_matrix)){
    P = PQ_com[ind_matrix[i, ], ]; Q = PQ_com[-ind_matrix[i, ], ]
    accuracy = c(accuracy, sum(c(sapply(classification[1: ntst], function(x){
      v = c(ip, iq, j) %in% x
      return(sum(colSums(t(P) == v) == ncol(P)))
    }), sapply(classification[-c(1: ntst)], function(x){
      v = c(ip, iq, j) %in% x
      return(sum(colSums(t(Q) == v) == ncol(Q)))
    })))/length(classification))
  }
  method_ind = which.max(accuracy)
  return(list("omega " = omega,"best_accuracy" = accuracy[method_ind]))
}
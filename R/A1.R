A1 <- function(data,tst_data,indice){
  library(lpSolve);library(gtools)
  accuracy = c()
  #give out omega based on input indice
  ip = indice[1]; jp = indice[2]; iq = indice[3]; jq = indice[4]
  n = nrow(data)/2; ntst = nrow(tst_data)/2
  obj = c(1, rep(0, 4), 
          c(rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n)))
  conp = rbind(cbind(rep(1, n), rep(-1, n), rep(1, n), rep(0, n), rep(0, n)), 
               cbind(rep(0, n), rep(-1, n), rep(1, n), rep(0, n), rep(0, n)), 
               cbind(rep(0, n), rep(0, n), rep(-1, n), rep(1, n), rep(0, n)),
               cbind(rep(0, n), rep(0, n), rep(-1, n), rep(0, n), rep(1, n)))
  conq = rbind(cbind(rep(1, n), rep(0, n), rep(0, n), rep(-1, n), rep(1, n)), 
               cbind(rep(0, n), rep(0, n), rep(0, n), rep(-1, n), rep(1, n)), 
               cbind(rep(0, n), rep(1, n), rep(0, n), rep(0, n), rep(-1, n)), 
               cbind(rep(0, n), rep(0, n), rep(1, n), rep(0, n), rep(-1, n)))
  con = cbind(rbind(conp, conq), diag(-1, nrow = 8*n, ncol = 8*n))
  rhs = c(rep(data[1: n, ip] - data[1: n, jp], 2), 
          data[1: n, jp] - data[1: n, iq], data[1: n, jp] - data[1: n, jq], 
          rep(data[-c(1: n), iq] - data[-c(1: n), jq], 2), 
          data[-c(1: n), jq] - data[-c(1: n), ip], 
          data[-c(1: n), jq] - data[-c(1: n), jp])
  dir=rep("<=",nrow(con))
  omega = rep(0, ncol(data))
  sol = lp("max", obj, con, dir, rhs)$solution
  omega[c(ip, jp, iq, jq)] = sol[2: 5]
  #if the dim of data is 5, there will be an error because dim(data[, -c(ip, jp, iq, jq)]) == null
  omega[-c(ip, jp, iq, jq)] = apply(data[, -c(ip, jp, iq, jq)], 
                                    2, 
                                    function(x){min(c(min(data[1: n, jp] - x[1:n])+omega[jp], 
                                                      min(data[-c(1: n), jq] - x[-c(1:n)])+omega[jq]))})
  classification = lapply(as.list(as.data.frame(apply(tst_data, 1, function(x){x+omega}))), 
                          function(x){which(x == max(x))})
  #choose classification standard
  P_base = matrix(c(1, 0, 0, 0,
                    0, 1, 0, 0,
                    1, 1, 0, 0, 
                    1, 1, 1, 1), ncol = 4, byrow = T); 
  Q_base = matrix(c(0, 0, 1, 0,
                    0, 0, 0, 1,
                    0, 0, 1, 1,
                    0, 0, 0, 0), ncol = 4, byrow = T); 
  PQ_com = matrix(c(1, 0, 1, 0,
                    1, 0, 0, 1,
                    0, 1, 1, 0,
                    0, 1, 0, 1,
                    1, 1, 1, 0,
                    1, 1, 0, 1,
                    1, 0, 1, 1,
                    0, 1, 1, 1), ncol = 4, byrow = T)
  ind_matrix = combinations(8, 4)
  for (l in 1:nrow(ind_matrix)){
    P = rbind(P_base, PQ_com[ind_matrix[l, ], ]); Q = rbind(Q_base, PQ_com[-ind_matrix[l, ], ])
    accuracy = c(accuracy, sum(c(sapply(classification[1: ntst], function(x){
      v = c(ip, jp, iq, jq) %in% x;
      return(sum(colSums(t(P) == v) == ncol(P)))
    }), sapply(classification[-c(1: ntst)], function(x){
      v = c(ip, jp, iq, jq) %in% x;
      return(sum(colSums(t(Q) == v) == ncol(Q)))
    })))/length(classification))
  }
  method_ind = which.max(accuracy)
  return(list("omega " = omega,"best_accuracy" = accuracy[method_ind]))
}



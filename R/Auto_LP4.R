ALP4 = function(data,CV = F){
  library(plyr);library(gtools);library(lpSolve)
  datasize = nrow(data)/2
  d = ncol(data)
  
  #cross validation
  CVgroup <- function(k,datasize){
    cvlist <- list()
    n <- rep(1:k,ceiling(datasize/k))[1:datasize]
    temp <- sample(n,datasize)
    x <- 1:k
    dataseq <- 1:datasize
    cvlist <- lapply(x,function(x) dataseq[temp==x])
    return(cvlist)
  }
  
  #5 fold CV
  k = 5
  cvlist1 = CVgroup(k,datasize) #P
  cvlist2 = CVgroup(k,datasize) #Q
  
  indice_box = permutations(d,3)
  if (CV == F){
    tst_data = rbind(data[cvlist1[[1]],],data[datasize+cvlist2[[1]],])
    data = data[-c(cvlist1[[1]],datasize+cvlist2[[1]]),]
    accu = lapply(1:nrow(indice_box),function(x){A4(data,tst_data,indice_box[x,])$best_accuracy})
    indice_num = which.max(accu)
    op_omega = A4(data,tst_data,indice_box[indice_num,])$omega
    op_accu = A4(data,tst_data,indice_box[indice_num,])$best_accuracy
    return(list("omega" = op_omega,"best_accuracy" = op_accu))
  }else{
    op_accu = 0
    op_omega = c()
    for (i in 1:k){
      tst_data = rbind(data[cvlist1[[i]],],data[(datasize+cvlist2[[i]]),])
      trn_data = data[-c(cvlist1[[i]],datasize+cvlist2[[i]]),]
      accu = lapply(1:nrow(indice_box),function(x){A4(trn_data,tst_data,indice_box[x,])$best_accuracy})
      indice_num = which.max(accu)
      omega = A4(trn_data,tst_data,indice_box[indice_num,])$omega
      accu = A4(trn_data,tst_data,indice_box[indice_num,])$best_accuracy
      if (accu > op_accu){
        op_accu = accu
        op_omega = omega
      }
    }
    return(list("omega" = op_omega,"best_accuracy" = op_accu))
  }
}

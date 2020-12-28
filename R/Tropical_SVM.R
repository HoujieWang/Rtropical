TSVM = function(data,CV = F,method = 1){
  methodbox = c(ALP1,ALP2,ALP3,ALP4)
  return(methodbox[[method]](data,CV))
}
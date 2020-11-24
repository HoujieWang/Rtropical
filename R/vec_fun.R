#' Vectorize a tree
#'
#' @importFrom ape cophenetic.phylo
#' @param x A tree object in R as list
#' @export
#' @examples
#' vec_fun(tree)
#'
vec_fun<-function(x){
  x = cophenetic.phylo(x)
  x[lower.tri(x)]
  # m<-dim(x)[1]
  # vecTreesVec<-rep(NA,choose(m,2))
  # for(row.num in 1:(m-1)){
  #   for(col.num in (row.num+1):m){
  #     vecTreesVec[col.num-row.num+(m-1+(m-1-row.num+2))*(row.num-1)/2]<-x[row.num,col.num]
  #   }
  # }
  # vecTreesVec
}

#' Compute Tropical Line Segment
#'
#' @param D1 a vector of length 2.
#' @param D2 a vector of length 2.
#' @param flag a logical value indicating if to swap \code{D1} and \code{D2}. (default: FALSE)
#'
#' @return A list containing the following:
#' \item{1}{End points as columns of each line component.}
#' \item{2}{Overall length of the tropical line segment}
#'
#' @author Qiwen Kang and Houjie Wang
#' Maintainer: Houjie Wang \email{wanghoujie6688@@gmail.com}
#'
#' @references Lin, B., Monod, A. and Yoshida, R. (2020)
#' \emph{Tropical Geometric Variation of Phylogenetic Tree Shapes}
#' \url{https://arxiv.org/pdf/2010.06158.pdf}
#'
#' @examples
#' seg <- tropseg(1: 2, 3: 4)
#' plot(x = seg[[1]][1, ], y = seg[[1]][2, ], "l")
#' points(x = seg[[1]][1, ], y = seg[[1]][2, ])
#' @export
#' @export tropseg
tropseg <- function(D1, D2, flag = FALSE){
  k <- length(D1)
  if(k != 2) warning("dimension has to be 2!")
  for(i in 1:k)
    D1[i] <- round(D1[i], 4)
  for(i in 1:k)
    D2[i] <- round(D2[i], 4)
  if(length(D2) != k)
    warning("dimension is wrong!")
  addd <- 0
  if(flag){
    tmp.D <- D2
    D2 <- D1
    D1 <- tmp.D
  }
  tmp.metric <- (D2 - D1)
  sorted.tmp.metric <- sort.int(tmp.metric, index.return=TRUE)
  D <- rep(0, k)

  D[sorted.tmp.metric$ix[2]] <- D2[sorted.tmp.metric$ix[2]]
  D[sorted.tmp.metric$ix[1]] <- min(D2[sorted.tmp.metric$ix[2]] - D1[sorted.tmp.metric$ix[2]] + D1[sorted.tmp.metric$ix[1]], D1[sorted.tmp.metric$ix[1]])

  distance <- max(abs(D1 - D))
  distance <- distance + max(abs(D2 - D))

  segment <- matrix(rep(0, 6), nrow=2, ncol=3)
  segment[,1] <- D1
  segment[,2] <- D
  segment[,3] <- D2

  return(list(segment, distance))
}

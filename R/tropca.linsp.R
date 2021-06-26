#' Tropical Principal Component Analysis by Tropical Linear Space
#'
#' Approximate the principal component as a tropical linear space
#' for a given data matrix and returns the results as an object of class \code{tropca}.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#'
#' @param x a data matrix, of size n x e, with each row an observation vector.
#' e is the dimension of the tropical space
#' @param pcs a numeric value indicating the order of principal component. (default: 2)
#' @param iteration a list with arguments controlling the iteration of the algorithm.
#' \describe{
#' \item{exhaust}{a logical variable indicating if to iterate over all possible combinations of the linear space
#' based on the given data matrix \code{x}. If FALSE, please input a number of iteration for \code{niter}.
#' If TRUE, please enter 0 for \code{niter} and this function will iterate over all possible combinations of linear space.
#' This could be time consuming when \code{x} is large. (default: FALSE)}
#' \item{niter}{a numeric variable indicating the number of iterations. (default: 100)}
#' }
#' @param ncores a numeric value indicating the number of threads utilized for multi-cored CPUs. (default: 2)
#'
#' @return A list of S3 class \code{"tropca"}, including:
#' \item{pc}{The principal component as a tropical linear space}
#' \item{obj}{The tropical PCA objective, the sum of tropical distance from each point to the projection.}
#' \item{projection}{The projections of all data points.}
#' \item{type}{The geometry of principal component.}
#'
#' @author Houjie Wang
#'
#' @references Yoshida, R.,Zhang, L. and Zhang, X. (2019)
#' \emph{Tropical Principal Component Analysis and Its Application to Phylogenetics,
#'  Bulletin of Mathematical Biology, 81, 568–597.}
#' \url{https://arxiv.org/pdf/2010.06158.pdf}.
#'
#' @keywords Tropical Geometry, Supervised Learning, Non-Euclidean Data
#' @examples
#' \dontrun{
#' library(Rfast)
#' n <- 100
#' e <- 10
#' sig2 <- 1
#' x <- rbind(rmvnorm(n, mu = c(5, -5, rep(0, e - 2)), sigma = diag(sig2, e)))
#' tropca_fit <- tropca.linsp(x)
#' }
#'
#' @export
#' @export tropca.linsp
tropca.linsp <- function(x, pcs = 2, iteration = list(), ncores = 2) {
  con <- list(
    exhaust = FALSE,
    niter = 100
  )
  con[names(iteration)] <- iteration

  exhaust <- con$exhaust
  niter <- con$niter
  pcs <- pcs + 1
  all_choices <- comboGeneral(nrow(x), pcs)
  if (exhaust) {
    all_choices <- lapply(1:nrow(all_choices), function(i) all_choices[i, ])
  } else {
    all_choices <- lapply(sample(1:nrow(all_choices), niter, replace = F), function(i) all_choices[i, ])
  }
  cl <- makeCluster(ncores)
  all_objs <- unlist(parLapply(cl, all_choices, function(ind) {
    V <- x[ind, ]
    proj <- troproj.linsp(x, V)
    temp <- x - proj
    sum(rowMaxs(temp, T) - rowMins(temp, T))
  }))
  stopCluster(cl)
  best_choice <- all_choices[[which.min(all_objs)]]
  pc <- x[best_choice, ]
  rownames(pc) <- paste("pc", 1:pcs, sep = "")
  proj_points <- troproj.linsp(x, pc)
  tropca.out <- list(
    "pc" = pc,
    "obj" = min(all_objs),
    "projection" = proj_points,
    "type" = "linear space"
  )
  class(tropca.out) <- "tropca"
  tropca.out
}

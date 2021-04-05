#' Tropical Principle Components Analysis by Tropical Linear Space
#'
#' Approximate the principle components as a tropical linear space for a given data matrix in the settings of
#' tropical geometry and return the results as an object of class tropca.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param pcs a numeric value indicating the number of principle components. (default: 3)
#' @param size a numeric value indicating the number of tropical linear spaces for where the principle components are searched. (default: 100)
#' @param ncores a numeric value indicating the number of threads utilized for multi-cored CPUs. (default: 2)
#' @param exhaust a logical variable indicating if to search
#'
#' @return A list of S3 class \code{"tropca"} containing the fitted model, including:
#' \item{pc}{The approximated principle components as a tropical polytope.}
#' \item{obj}{The tropical PCA objective, the sum of tropical distances from points to their projections on the principle components.}
#' \item{projection}{The projections of points.}
#'
#' @author Houjie Wang
#'
#' @references Page, R., Yoshida, R. & Zhang L.
#' \emph{Tropical principal component analysis on the space of phylogenetic trees.
#' J. Bioinform., Volume 36, Issue 17, 4590–4598 (2020).}
#' \url{https://doi.org/10.1093/bioinformatics/btaa564}
#'
#' @keywords Tropical Geometry, Supervised Learning, Non-Euclidean Data
#' @examples
#' \dontrun{
#' library(Rfast)
#' n <- 100; e <- 10; sig2 <- 1; pcs <- 3
#' x <- rbind(rmvnorm(n, mu = c(5, -5, rep(0, e-2)), sigma = diag(sig2, e)))
#' tropca_fit <- tropca.linsp(x)
#' }
#'
#' @export
#' @export tropca.linsp
tropca.linsp = function(x, pcs = 3, size = 100, ncores = 2, exhaust = FALSE){
  all_choices <- comboGeneral(nrow(x), pcs)
  if (exhaust){
    all_choices <- lapply(1: nrow(all_choices), function(i) all_choices[i, ])
  } else{
    all_choices <- lapply(sample(1: nrow(all_choices), size, replace = F), function(i) all_choices[i, ])
  }
  cl <- makeCluster(ncores)
  all_objs <- unlist(parLapply(cl, all_choices, function(ind){
    V <- x[ind, ]
    proj <- troproj.linsp(x, V)
    temp <- x - proj
    sum(rowMaxs(temp, T) - rowMins(temp, T))
  }))
  best_choice <- all_choices[[which.min(all_objs)]]
  pc <- x[best_choice, ]
  rownames(pc) <- paste("pc", 1: pcs, sep = "")
  proj_points <- troproj.linsp(x, pc)
  tropca.out <- list("pc" = pc,
                     "obj" = min(all_objs),
                     "projection" = proj_points)
  class(tropca.out) <- "tropca"
  tropca.out
}

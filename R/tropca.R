#' Tropical Principle Components Analysis
#'
#' Approximate the principle components as a polytope for a given data matrix in the settings of
#' tropical geometry via MCMC and return the results as an object of class tropca.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel makeCluster
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#' @importFrom stats runif
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param pcs a numeric value indicating the number of principle components. (default: 3)
#' @param nsample a numeric value indicating the number of samples of MCMC. (default: 1000)
#' @param ncores a numeric value indicating the number of threads utilized for multi-cored CPUs. (default: 2)
#'
#' @return A list of S3 class \code{"tropca"} containing the fitted model, including:
#' \item{pc}{The approximated principle components as a tropical polytope.}
#' \item{obj}{The tropical PCA objective, the sum of tropical distances from points to their projections on the principle components.}
#' \item{projection}{The projections of points.}
#' \item{samples}{All sampled tropical polytopes.}
#' \item{objs}{All tropical PCA objectives for all sampled tropical polytopes.}
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
#' n <- 50; e <- 50; s <- 5
#' x <- rbind(rmvnorm(n, mu = c(5, -5, rep(0, e-2)), sigma = diag(s, e)),
#'            rmvnorm(n, mu = c(-5, 5, rep(0, e-2)), sigma = diag(s, e)))
#' tropca_fit <- tropca(x)
#' plot(tropca_fit)
#' }
#'
#' @export
#' @export tropca

tropca <- function(x, pcs = 3, nsample = 1000, ncores = 2){
  n = nrow(x)
  cl <- makeCluster(ncores)
  x_list <- lapply(seq_len(n), function(i) x[i, ])
  tropca_objs = vector(mode = "numeric", nsample)
  samples = matrix(NA, nrow = nsample, ncol = pcs)
  samples[1, ] = sample(1: n, pcs)
  tropca_objs[1] = tropca.obj(t(x[samples[1, ], ]), x_list, cl)

  t = 1
  while (t < nsample){
    # Find a new proposal by changing a randomly selected vertex of the current polytope
    current_choice = samples[t, ]
    current_obj = tropca_objs[t]

    change_ind <- sample(pcs, 1)
    out_change <- sample(c(1: n)[-current_choice], 1)
    new_choice <- c(current_choice[-change_ind], out_change)
    new_obj = tropca.obj(t(x[new_choice, ]), x_list, cl)

    # Compute the probability we accept the new PCA base
    p = min(1, current_obj/new_obj)

    if(sample(c(0, 1), 1, prob = c(1 - p, p)) == 1){
      samples[(t + 1), ] <- new_choice
      tropca_objs[(t + 1)] <- new_obj
      t = t + 1
    }
  }
  min_index <- which(tropca_objs == min(tropca_objs))[1]
  best_obj = tropca_objs[min_index]
  pc = x[samples[min_index, ], ]
  proj_points <- do.call("rbind", parLapply(cl, x_list, troproj , tconv = t(pc)))
  rownames(pc) <- paste("pc", 1: pcs, sep = "")
  tropca.out <- list("pc" = pc,
                     "obj" = tropca_objs[min_index],
                     "projection" = proj_points,
                     "samples" = samples,
                     "objs" = tropca_objs)
  class(tropca.out) <- "tropca"
  tropca.out
}

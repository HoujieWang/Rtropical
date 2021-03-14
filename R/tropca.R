#' Tropical Principle Components Analysis
#'
#' Approximate the principle components as a polytope for a given data matrix in the settings of
#' tropical geometry via MCMC and returns the results as an object of class tropca.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel makeCluster
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#' @importFrom stats runif
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param pcs a numeric value indicating the number of principle components. (default: 3)
#' @param rep a numeric value indicating the number of repetitions of MCMC. (default: 2)
#'
#' @return An object with S3 class \code{"tropca"} containing the fitted model, including:
#' \item{pc}{The principle components.}
#' \item{obj}{The sum of tropical distances from each point to its projection on the principle components.}
#' \item{projection}{The projections of given data set.}
#'
#' @author Qiwen Kang and Houjie Wang
#'
#' Maintainer: Houjie Wang \email{wanghoujie6688@@gmail.com}
#' @references Page, R., Yoshida, R. & Zhang L.
#' \emph{Tropical principal component analysis on the space of phylogenetic trees.
#' J. Bioinform., Volume 36, Issue 17, 4590–4598 (2020).}
#' \url{https://doi.org/10.1093/bioinformatics/btaa564}
#'
#' @keywords Tropical Geometry, Supervised Learning, Non-Euclidean Data
#' @examples
#' library(Rfast)
#' n <- 20; e <- 50; s <- 5
#' x <- rbind(rmvnorm(n, mu = c(5, -5, rep(0, e-2)), sigma = diag(s, e)),
#'            rmvnorm(n, mu = c(-5, 5, rep(0, e-2)), sigma = diag(s, e)))
#' tropca_fit <- tropca(x)
#' plot(tropca_fit)
#' @export
#' @export tropca

tropca <- function(x, pcs = 3, rep = 2){
  tropca_objs <- rep(NA, rep)
  comb_list <- list()
  cl <- makeCluster(2)
  N <- nrow(x)
  x_list <- lapply(seq_len(nrow(x)), function(i) x[i, ])
  for(j in 1:rep){
    sample_init <- sample(N, pcs)
    best <- 100000
    out <- c(1:N)[-sample_init]
    pc_base_init <- x[sample_init, ]
    init_value <- tropca.obj(t(pc_base_init), x_list, cl)
    t = 0
    while(length(out)!=0){
      t = t+ 1
      change_ind <- sample(pcs,1)
      out_change <- sample(out, 1)
      comb_set <- c(sample_init[-change_ind], out_change)

      new_base <- x[comb_set, ]

      update_value <- tropca.obj(t(new_base), x_list, cl)
      r <- init_value/update_value

      if(runif(1) < min(r, 1)){
        sample_init <- comb_set

        best <- ifelse(update_value < best, update_value, best)
      }
      out <- out[-which(out==out_change)]
      init_value <- update_value
    }
    tropca_objs[j] <- best
    comb_list[[j]] <- sample_init
  }
  min_index <- which(tropca_objs==min(tropca_objs))
  pc <- x[comb_list[[min_index]], ]
  proj_points <- do.call("rbind", parLapply(cl, x_list, troproj , tconv = t(pc)))
  rownames(pc) <- paste("pc", 1: pcs, sep = "")
  tropca.out <- list("pc" = pc,
                     "obj" = tropca_objs[min_index],
                     "projection" = proj_points)
  class(tropca.out) <- "tropca"
  tropca.out
}

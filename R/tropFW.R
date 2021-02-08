#' Find the Tropical Fermat-Weber Point
#'
#' Funcition to compute a tropical fermar-weber point such that the sum of tropical distances is minimized.
#'
#' @importFrom RcppAlgos comboGeneral
#' @importFrom lpSolve lp
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#'
#' @return A list containing the optimal fermat-weber point and the minimized sum of distance:
#' \item{fw}{The optimal fermat-weber point.}
#' \item{distsum}{The sum of distance from each observation to the optimal fermat-weber point.}
#'
#' @author Houjie Wang
#'
#' Maintainer: Houjie Wang \email{wanghoujie6688@@gmail.com}
#'
#' @references Lin, B.,Sturmfels, B., Tang, X. and Yoshida, R. (2015)
#' \emph{Convexity in Tree Spaces, SIAM Journal on Discrete Mathematics, Vol. 31(3) 2015–2038}
#'
#' @keywords Tropical Geometry, Fermat-Weber Points
#'
#' @examples
#'
#' x = matrix(rnorm(100), ncol = 10)
#' tropFW(x)
#'
#' @export
#' @export tropFW
tropFW <- function(x){
  nn <- nrow(x)
  e <- ncol(x)
  jk <- comboGeneral(1: e, 2)
  combn_size <- nrow(jk)
  conY <- matrix(0, nrow = nn*combn_size, ncol = e)
  all_v <- x[, jk[, 2]] - x[, jk[, 1]]
  conY[cbind(1: (nn*combn_size), rep(jk[, 1], each = nrow(x)))] <- -1
  conY[cbind(1: (nn*combn_size), rep(jk[, 2], each = nrow(x)))] <- 1
  conD <- matrix(0, nrow = 2*nn*combn_size, ncol = nn)
  conD[cbind(1: (2*nn*combn_size), 1: nn)] <- -1
  con <- cbind(conD, rbind(conY, -conY))
  rhs <- c(matrix(all_v), -matrix(all_v))
  obj <- c(rep(1, nn), rep(0, e))
  dir <- rep("<=", (2*nn*combn_size))
  sol <- lp("min", obj, con, dir, rhs)
  list("fw" = sol$solution[-c(1: nn)], "distsum" = sol$objval)
}

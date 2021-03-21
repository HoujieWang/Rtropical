#' Tropical Fermat-Weber Point
#'
#' Compute the tropical fermat-weber point for a given data set such that the sum of tropical distance to each point is minimized.
#'
#' @importFrom RcppAlgos comboGeneral
#' @importFrom lpSolveAPI make.lp
#' @importFrom lpSolveAPI set.constr.type
#' @importFrom lpSolveAPI set.rhs
#' @importFrom lpSolveAPI solve.lpExtPtr
#' @importFrom lpSolveAPI get.variables
#' @importFrom lpSolveAPI set.column
#' @importFrom lpSolveAPI get.objective
#' @importFrom lpSolveAPI set.objfn
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#'
#' @return A list containing:
#' \item{fw}{The fermat-weber point.}
#' \item{distsum}{The sum of distance from each observation to the fermat-weber point.}
#'
#' @author Houjie Wang
#'
#' Maintainer: Houjie Wang \email{wanghoujie6688@@gmail.com}
#'
#' @references Lin, B.,Sturmfels, B., Tang, X. and Yoshida, R. (2015)
#' \emph{Convexity in Tree Spaces, SIAM Journal on Discrete Mathematics, Vol. 31(3) 2015–2038}
#'
#' @keywords Tropical Geometry, Fermat-Weber Point
#'
#' @examples
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
  obj <- c(rep(1, nn), rep(0, e))
  conY <- matrix(0, nrow = nn*combn_size, ncol = e)
  all_v <- x[, jk[, 2]] - x[, jk[, 1]]
  conY[cbind(1: (nn*combn_size), rep(jk[, 1], each = nrow(x)))] <- -1
  conY[cbind(1: (nn*combn_size), rep(jk[, 2], each = nrow(x)))] <- 1
  conD <- matrix(0, nrow = 2*nn*combn_size, ncol = nn)
  conD[cbind(1: (2*nn*combn_size), 1: nn)] <- -1
  con <- cbind(conD, rbind(conY, -conY))
  conD_2 <- matrix(0, nrow = nn*e, ncol = nn)
  conD_2[cbind(1: (nn*e), 1: nn)] <- -1
  conY_2 <- matrix(0, nrow = nn*e, ncol = e)
  conY_2[cbind(1: (nn*e), rep(1: e, each = nn))] <- 1
  con_2 <- cbind(conD_2, conY_2)
  con_new <- rbind(con, con_2)
  rhs_new <- c(matrix(all_v), -matrix(all_v), rep(0, (nn*e)))
  dir_new <- rep("<=", (2*nn*combn_size + (nn*e)))

  lprec <- make.lp(nrow(con_new), nn+e)
  for (i in 1: ncol(con_new)){set.column(lprec, i, -con_new[, i])}
  set.constr.type(lprec, rep(">=", (2*nn*combn_size + (nn*e))))
  set.rhs(lprec, -rhs_new)
  set.objfn(lprec, c(rep(1, nn), rep(0, e)))
  solve.lpExtPtr(lprec)
  sols = get.variables(lprec)
  list("fw" = sols[-c(1: nn)], "distsum" = sum(sols[1: nn]))
}

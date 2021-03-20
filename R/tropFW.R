#' Tropical Fermat-Weber Point
#'
#' Compute the tropical fermat-weber point for a given data set such that the sum of tropical distance to each point is minimized.
#'
#' @importFrom RcppAlgos comboGeneral
#' @importFrom lpSolveAPI make.lp
#' @importFrom lpSolveAPI add.constraint
#' @importFrom lpSolveAPI solve.lpExtPtr
#' @importFrom lpSolveAPI get.variables
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
  n <- dim(x)[1]
  m <- dim(x)[2]
  lprec <- make.lp(0, n+m)
  objective <- mat.or.vec(n+m,1)
  for (i in seq(n)) {
    objective[i] <- 1
  }
  set.objfn(lprec, objective)
  for (i in seq(n)) {
    for (j in seq(m)) {
      for (k in seq(m)) {
        v <- mat.or.vec(n+m,1)
        v[i] <- 1
        v[n+k] <- 1
        v[n+j] <- -1
        add.constraint(lprec, v, ">=", x[i,k] - x[i,j])
      }
    }
  }
  solve.lpExtPtr(lprec)
  sols <- get.variables(lprec)
  list("fw" = sols[-c(1: n)], "distsum" = get.objective(lprec))
}

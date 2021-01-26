#' Extract Optimal Tropical Hyperplane from a cv.tropsvm object
#'
#' Obtain the optimal tropical hyperplane in the form of vectors from a cv.tropsvm object.
#'
#' @param object a fitted \code{"cv.tropsvm"} object.
#'
#' @return An output of the apex of the fitted optimal tropical hyperplane.
#' @method coef cv.tropsvm
#' @export
#' @export coef.cv.tropsvm
coef.cv.tropsvm <- function(object, ...){
  cat("The apex of seprarating hyperplane is: \n", -object$coef)
}

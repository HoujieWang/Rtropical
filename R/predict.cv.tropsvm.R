#' Predict Method for Tropical Support Vector Machines based on Cross-Validation
#'
#' Predicts values based upon a model trained by \code{cv.tropsvm}.
#' @importFrom RcppAlgos comboGeneral
#'
#' @param object a fitted \code{"cv.tropsvm"} object.
#' @param newx a data matrix, of dimension nobs x nvars used as testing data.
#' @param newy a response vector with one label for each row/component of x used as testing data.
#'
#' @return A vector of predicted values of a vector of labels.
#'
#' @author Houjie Wang
#' Maintainer: Houjie Wang \email{whj666@@uw.edu}
#'
#' @references Tang, X., Wang, H. and Yoshida, R. (2020)
#' \emph{Tropical Support Vector Machine and its Applications to Phylogenomics}
#' \url{https://arxiv.org/pdf/2003.00677.pdf}
#'
#' @seealso \code{summary}, \code{coef} and the \code{cv.tropsvm} function.
#'
#' @examples
#'
#' # data generation
#' library(Rfast)
#' e <- 100; n = 100; N = 100; s = 10
#' x <- rbind(rmvnorm(n, mu = c(5, -5, rep(0, e-2)), sigma = diag(s, e)),
#'           rmvnorm(n, mu = c(-5, 5, rep(0, e-2)), sigma = diag(s, e)))
#' y <- as.factor(c(rep(1, n), rep(2, n)))
#' newx <- rbind(rmvnorm(N, mu = c(5, -5, rep(0, e-2)), sigma = diag(s, e)),
#'              rmvnorm(N, mu = c(-5, 5, rep(0, e-2)), sigma = diag(s, e)))
#' newy <- as.factor(rep(c(1, 2), each = N))
#'
#' # train the tropical svm
#' cv_tropsvm_fit <- cv.tropsvm(x, y)
#'
#' # test with new data
#' pred <- predict(cv_tropsvm_fit , newx, newy)
#'
#' # check with accuracy
#' table(pred, newy)
#'
#' # compute testing accuracy
#' sum(predict(cv_tropsvm_fit , newx, newy) == newy)/length(newy)
#' @method predict cv.tropsvm
#'
#' @export
#' @export predict.cv.tropsvm
predict.cv.tropsvm <- function(object, newx, newy, ...){
  P_base <- matrix(c(1, 0, 0, 0,
                     0, 1, 0, 0,
                     1, 1, 0, 0,
                     1, 1, 1, 1), ncol = 4, byrow = T);
  Q_base <- matrix(c(0, 0, 1, 0,
                     0, 0, 0, 1,
                     0, 0, 1, 1,
                     0, 0, 0, 0), ncol = 4, byrow = T);
  PQ_com <- matrix(c(1, 0, 1, 0,
                     1, 0, 0, 1,
                     0, 1, 1, 0,
                     0, 1, 0, 1,
                     1, 1, 1, 0,
                     1, 1, 0, 1,
                     1, 0, 1, 1,
                     0, 1, 1, 1), ncol = 4, byrow = T)
  colnames(PQ_com) <- c("ip", "jp", "iq", "jq")
  all_method_ind <- RcppAlgos::comboGeneral(8, 4)
  if (length(unique(newy)) != 2){
    stop("Only two classes are allowded.")
  }
  if (is.data.frame(x)){
    x <- data.matrix(x)
  }
  classes <- unique(newy)
  best_method <- object$`method index`
  omega <- object$coef
  best_assignment <- object$assignment[c(1, 3, 2, 4)]

  classification_method <- rbind(rbind(P_base, PQ_com[all_method_ind[best_method, ], ]),
                                 rbind(Q_base, PQ_com[-all_method_ind[best_method, ], ]))
  shifted_tst_data <- eachrow(newx, omega, "+")
  diff <- eachrow(t(shifted_tst_data), rowMaxs(shifted_tst_data, T), oper = "-")
  classification <- lapply(lapply(seq_len(ncol(diff)), function(i) diff[, i]), function(x){which(abs(x) < 1e-10)})
  classification <- sapply(classification, function(x){which(colSums(abs(t(classification_method) - best_assignment %in% x)) == 0)})
  classification[classification <= 8] <- classes[1]
  classification[classification > 8] <- classes[2]
  as.factor(classification)
}

#' Predict Method for Tropical Support Vector Machines
#'
#' Predicts values based upon a model trained by \code{tropsvm}.
#'
#' @importFrom RcppAlgos comboGeneral
#'
#' @param object a fitted \code{"cv.tropsvm"} object.
#' @param newx a data matrix, of dimension nobs x nvars used as testing data.
#' @param \dots Not used. Other arguments to predict.
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
#' @seealso \code{summary}, \code{coef} and the \code{tropsvm} function.
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
#' tropsvm_fit <- tropsvm(x, y, auto.assignment = TRUE, ind = 1)
#'
#' # test with new data
#' pred <- predict(tropsvm_fit , newx)
#'
#' # check with accuracy
#' table(pred, newy)
#'
#' # compute testing accuracy
#' sum(pred == newy)/length(newy)
#' @method predict tropsvm
#' @export
#' @export predict.tropsvm
predict.tropsvm <- function(object, newx, ...){
  # object = tropsvm_fit;
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
  if (is.data.frame(x)){
    x <- data.matrix(x)
  }
  classes <- object$`levels`
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

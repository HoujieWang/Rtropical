#' Predict Method for Tropical Support Vector Machines based on Cross-Validation
#'
#' Predicts values based upon a model trained by \code{cv.tropsvm}.
#' @importFrom RcppAlgos comboGeneral
#'
#' @param object a fitted \code{"cv.tropsvm"} object.
#' @param newx a data matrix, of dimension nobs x nvars used as testing data.
#' @param \dots Not used. Other arguments to predict.

#' @return A vector of predicted values of a vector of labels.
#'
#'
#' @seealso \code{summary}, \code{coef} and the \code{cv.tropsvm} function.
#'
#' @examples
#'
#' # data generation
#' library(Rfast)
#' e <- 20
#' n <- 10
#' N <- 10
#' s <- 5
#' x <- rbind(
#'   rmvnorm(n, mu = c(5, -5, rep(0, e - 2)), sigma = diag(s, e)),
#'   rmvnorm(n, mu = c(-5, 5, rep(0, e - 2)), sigma = diag(s, e))
#' )
#' y <- as.factor(c(rep(1, n), rep(2, n)))
#' newx <- rbind(
#'   rmvnorm(N, mu = c(5, -5, rep(0, e - 2)), sigma = diag(s, e)),
#'   rmvnorm(N, mu = c(-5, 5, rep(0, e - 2)), sigma = diag(s, e))
#' )
#' newy <- as.factor(rep(c(1, 2), each = N))
#'
#' # train the tropical svm
#' cv_tropsvm_fit <- cv.tropsvm(x, y, parallel = FALSE)
#'
#' # test with new data
#' pred <- predict(cv_tropsvm_fit, newx)
#'
#' # check with accuracy
#' table(pred, newy)
#'
#' # compute testing accuracy
#' sum(pred == newy) / length(newy)
#' @method predict cv.tropsvm
#'
#' @export
#' @export predict.cv.tropsvm
predict.cv.tropsvm <- function(object, newx, ...) {
  # object = cv.svmmodel;
  classes <- object$`levels`
  best_method <- object$`index`
  omega <- object$apex
  best_assignment <- object$assignment[c(1, 3, 2, 4)]
  ip <- best_assignment[1]
  jp <- best_assignment[2]
  iq <- best_assignment[3]
  jq <- best_assignment[4]

  if (length(unique(best_assignment)) == 2) {
    shifted_tst_data <- eachrow(newx, omega, "+")
    classification <- rowMaxs(shifted_tst_data)
    classification[classification == ip] <- levels(classes)[1]
    classification[classification == iq] <- levels(classes)[2]
    as.factor(classification)
  } else {
    all_method_ind <- comboGeneral(8, 4)
    if (length(unique(best_assignment)) == 4) {
      P_base <- matrix(c(
        1, 0, 0, 0,
        0, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 1, 1
      ), ncol = 4, byrow = T)
      Q_base <- matrix(c(
        0, 0, 1, 0,
        0, 0, 0, 1,
        0, 0, 1, 1,
        0, 0, 0, 0
      ), ncol = 4, byrow = T)
      PQ_com <- matrix(c(
        1, 0, 1, 0,
        1, 0, 0, 1,
        0, 1, 1, 0,
        0, 1, 0, 1,
        1, 1, 1, 0,
        1, 1, 0, 1,
        1, 0, 1, 1,
        0, 1, 1, 1
      ), ncol = 4, byrow = T)
    }
    if (length(unique(best_assignment)) == 3) {
      P_base <- c()
      Q_base <- c()
      PQ_com <- matrix(c(
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        1, 1, 0,
        1, 0, 1,
        0, 1, 1,
        1, 1, 1,
        0, 0, 0
      ), ncol = 3, byrow = T)
      if (ip == jq) {
        PQ_com <- PQ_com[, c(1, 2, 3, 1)]
      }
      if (iq == jp) {
        PQ_com <- PQ_com[, c(1, 2, 2, 3)]
      }
      if (jp == jq) {
        PQ_com <- PQ_com[, c(1, 2, 3, 2)]
      }
    }
    colnames(PQ_com) <- c("ip", "jp", "iq", "jq")
    classification_method <- rbind(
      rbind(P_base, PQ_com[all_method_ind[best_method, ], ]),
      rbind(Q_base, PQ_com[-all_method_ind[best_method, ], ])
    )
    shifted_tst_data <- eachrow(newx, omega, "+")
    diff <- eachrow(t(shifted_tst_data), rowMaxs(shifted_tst_data, T), oper = "-")
    classification <- lapply(lapply(seq_len(ncol(diff)), function(i) diff[, i]), function(x) {
      which(abs(x) < 1e-10)
    })
    classification <- sapply(classification, function(x) {
      which(colSums(abs(t(classification_method) - best_assignment %in% x)) == 0)
    })
    classification_temp = classification
    classification[classification_temp <= nrow(classification_method) / 2] <- levels(classes)[1]
    classification[classification_temp > nrow(classification_method) / 2] <- levels(classes)[2]
    names(classification) = NULL
    as.factor(classification)
  }
}

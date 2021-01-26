#' Summarize an Analysis of Cross-Validated Tropical Support Vector Machine
#'
#' Return a summary with a more detailed explanation of the object \code{"cv.tropsvm"}.
#'
#' @importFrom RcppAlgos comboGeneral
#'
#' @param object a fitted \code{"cv.tropsvm"} object
#'
#' @return A summary of the crucial information of a tropical support vector machine is printed, incluiding
#' the selected best assignment and classification methods and the validation accuracy of each data fold. The
#' summary section of classification methods specifies the sectors and their intersections used to classify
#' points of two different categories.
#'
#' @author Houjie Wang
#' Maintainer: Houjie Wang \email{whj666@@uw.edu}
#'
#' @references Tang, X., Wang, H. and Yoshida, R. (2020)
#' \emph{Tropical Support Vector Machine and its Applications to Phylogenomics}
#' \url{https://arxiv.org/pdf/2003.00677.pdf}
#'
#' @seealso \code{predict}, \code{coef} and the \code{cv.tropsvm} function.
#'
#' @examples
#'
#' # data generation
#' library(Rfast)
#' e <- 100; n = 100; N = 100; s = 10
#' x <- rbind(rmvnorm(n, mu = c(5, -5, 3, 10, rep(0, e-4)), sigma = diag(s, e)),
#'           rmvnorm(n, mu = c(-5, 5, -3, -10, rep(0, e-4)), sigma = diag(s, e)))
#' y <- as.factor(c(rep(1, n), rep(2, n)))
#' newx <- rbind(rmvnorm(N, mu = c(5, -5, rep(0, e-2)), sigma = diag(s, e)),
#'              rmvnorm(N, mu = c(-5, 5, rep(0, e-2)), sigma = diag(s, e)))
#' newy <- as.factor(rep(c(1, 2), each = N))
#'
#' # train the tropical svm with cross-validation
#' cv_tropsvm_fit <- cv.tropsvm(x, y)
#'
#' summary(cv_tropsvm_fit)
#' @method summary cv.tropsvm
#'
#' @export
#' @export summary.cv.tropsvm
summary.cv.tropsvm <- function(object, ...){
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
  best_method <- object$`method index`
  best_assignment <- object$`assignment`
  classification_method <- list("P method" = rbind(P_base, PQ_com[all_method_ind[best_method, ], ]), "Q method" = rbind(Q_base, PQ_com[-all_method_ind[best_method, ], ]))
  cat("Tropical SVM under ", object$nfold, "-fold cross validation: \n\n\n", sep = "")
  cat("Best assignment: ", paste(c("ip =", "jp =", "iq =", "jq ="), best_assignment, collapse = ", "), ".\n\n", sep = "")
  cat("Best classification method: \n\n")
  for (i in 1: nrow(classification_method[[1]])){
    row_i <- classification_method[[1]][i, ]
    if(sum(row_i) == 1){
      cat(best_assignment[which(row_i != 0)], "-th sector.\n", sep = "")
    } else{
      cat("The common boundary of", best_assignment[which(row_i != 0)], "-th sectors.\n")
    }
  }
  cat("Points on the locations above will be classified as : ", object$levels[1], ".\n\n", sep = "")

  for (i in 1: nrow(classification_method[[2]])){
    row_i <- classification_method[[2]][i, ]
    if(sum(row_i) == 1){
      cat(best_assignment[which(row_i != 0)], "-th sector.\n", sep = "")
    }
    if (sum(row_i) > 1){
      cat("The common boundary of", best_assignment[which(row_i != 0)], "-th sectors.\n")
    }
  }
  cat("Points on the locations above will be classified as : ", object$levels[2], ".\n\n", sep = "")
  cat("Best validation accuracy of each fold: ", paste(round(object$accuracy*100, digits = 4), "%", sep = ""), ".\n")
}

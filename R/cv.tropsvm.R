#' Cross-Validation for Tropical Support Vector Machines
#'
#' Conduct k-fold cross validation for tropsvm and return an object \code{"cv.tropsvm"}.
#'
#' @importFrom parallel parLapply
#' @importFrom parallel makeCluster
#' @importFrom parallel setDefaultCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom Rfast eachrow
#' @importFrom Rfast rowMaxs
#' @importFrom Rfast colMins
#' @importFrom lpSolve lp
#' @importFrom caret createFolds
#'
#' @param x a data matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y a response vector with one label for each row/component of x.
#' @param parallel a logical value indicating if parallel computing should be used. (default: TRUE)
#' @param nfold a numeric value of the number of data folds for cross-validation. (default: 10)
#' @param nassignment a numeric value indicating the size of the parameter grid of assignments. (default: 10)
#' @param ncores a numeric value indicating the number of threads utilized for multi-cored CPUs. (default: 2)
#'
#' @return An object with S3 class \code{"cv.tropsvm"}.
#' \item{coef}{The vector of the fitted optimal tropical hyperplane.}
#' \item{assignment}{The best \code{assignment} tuned by cross-validation.}
#' \item{method index}{The best classification method indices \code{ind} tuned by cross-validation.}
#' \item{levels}{The name of each category, consistent with categories in \code{y}.}
#' \item{accuracy}{The validation accuracy for each validation set.}
#' \item{nfold}{The number of folds used in cross-validation}
#'
#' @author Houjie Wang and Kaizhang Wang
#' Maintainer: Houjie Wang \email{whj666@@uw.edu}
#'
#' @references Tang, X., Wang, H. and Yoshida, R. (2020)
#' \emph{Tropical Support Vector Machine and its Applications to Phylogenomics}
#' \url{https://arxiv.org/pdf/2003.00677.pdf}
#'
#' @seealso \code{summary}, \code{predict}, \code{coef} and the \code{tropsvm} function.
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
#' summary(cv_tropsvm_fit)
#' coef(cv_tropsvm_fit)
#'
#' # test with new data
#' pred <- predict(cv_tropsvm_fit, newx)
#'
#' # check with accuracy
#' table(pred, newy)
#'
#' # compute testing accuracy
#' sum(pred == newy)/length(newy)
#'
#' @export
#' @export cv.tropsvm

cv.tropsvm <- function(x, y, parallel = TRUE, nfold = 10, nassignment = 10, ncores = 2){
  if (length(unique(y)) != 2){
    stop("Only two classes are allowded.")
  }
  if (is.data.frame(x)){
    x <- data.matrix(x)
    warning("Input data has to be of class 'Matrix'.")
  }

  classes <- unique(y)
  reorder_ind <- c(which(y == classes[1]), which(y == classes[2]))
  np = sum(y == classes[1]); nq = sum(y == classes[1])
  y <- y[reorder_ind]
  x <- x[reorder_ind, ]


  train_index <- caret::createFolds(y, k = nfold, returnTrain = TRUE)
  all_assignment <- matrix(0, nrow = nfold*nassignment, ncol = 4)
  for (i in 1: length(train_index)){
    P = x[train_index[[i]] <= np, ]
    Q = x[train_index[[i]] > np, ]
    all_assignment[((i-1)*nassignment+1): (i*nassignment), ] = assignment_finder(P, Q)[1: nassignment, ]
  }
  all_assignment = unique(all_assignment)
  all_assignment_list = lapply(seq_len(nrow(all_assignment)), function(i) all_assignment[i, ])
  all_accuracy_list <- list()
  for (i in 1: length(train_index)){
    # i = 1
    data <- x[train_index[[i]], ]
    label <- y[train_index[[i]]]
    n1 <- sum(label == classes[1])
    n2 <- sum(label == classes[2])
    n <- n1 + n2
    val_data <- x[-train_index[[i]], ]
    val_label <- y[-train_index[[i]]]
    val_n1 <- sum(val_label == classes[1])
    val_n2 <- sum(val_label == classes[2])
    val_n <- val_n1 + val_n2


    f.obj <- c(1, rep(0, 4), c(rep(-1, n1), rep(-1, n1), rep(-1, n1), rep(-1, n1), rep(-1, n2),
                               rep(-1, n2), rep(-1, n2), rep(-1, n2)))
    f.conp <- rbind(cbind(rep(1, n1), rep(-1, n1), rep(1, n1), rep(0, n1), rep(0, n1)),
                    cbind(rep(0, n1), rep(-1, n1), rep(1, n1), rep(0, n1), rep(0, n1)),
                    cbind(rep(0, n1), rep(0, n1), rep(-1, n1), rep(1, n1), rep(0, n1)),
                    cbind(rep(0, n1), rep(0, n1), rep(-1, n1), rep(0, n1), rep(1, n1)))
    f.conq <- rbind(cbind(rep(1, n2), rep(0, n2), rep(0, n2), rep(-1, n2), rep(1, n2)),
                    cbind(rep(0, n2), rep(0, n2), rep(0, n2), rep(-1, n2), rep(1, n2)),
                    cbind(rep(0, n2), rep(1, n2), rep(0, n2), rep(0, n2), rep(-1, n2)),
                    cbind(rep(0, n2), rep(0, n2), rep(1, n2), rep(0, n2), rep(-1, n2)))
    f.con <- cbind(rbind(f.conp, f.conq), diag(-1, nrow = 4*n, ncol = 4*n))
    f.dir <- rep("<=", n)

    if (parallel){
      cl <- makeCluster(ncores)
      clusterExport(cl, list("data", "label", "val_data", "val_label", "tropsvm", "lp", "eachrow", "rowMaxs", "colMins"), envir = environment())
      all_accuracy <- parLapply(cl, all_assignment_list, function(assignment){
        tropsvm(x = data, y = label, accuracy = T, assignment = 1: 4, ind = 1: 70, newx = val_data, newy = val_label)})
      stopCluster(cl)
    } else{
      all_accuracy <- lapply(all_assignment_list, function(assignment){tropsvm(data, label, accuracy = T, assignment, ind = 1: 70, newx = val_data, newy = val_label)})
    }
    accuracy_mat <- do.call("rbind", all_accuracy)
    all_accuracy_list[[i]] <- accuracy_mat
  }
  all_accuracy_mat <- Reduce("+", all_accuracy_list)
  best_hyperparms <- matrix(which(all_accuracy_mat == max(all_accuracy_mat), arr.ind = T)[1, ], ncol = 2, byrow = TRUE)
  best_assignment <- all_assignment[best_hyperparms[1, 1], ]
  best_method_ind <- best_hyperparms[1, 2]
  best_accuracy <- sapply(all_accuracy_list, max)
  best_fold <- which.max(sapply(all_accuracy_list, function(x){x[best_hyperparms]}))
  data <- x[train_index[[best_fold]], ]
  label <- y[train_index[[best_fold]]]
  n1 <- sum(label == classes[1])
  n2 <- sum(label == classes[2])
  n <- n1 + n2
  val_data <- x[-train_index[[best_fold]], ]
  val_label <- y[-train_index[[best_fold]]]
  val_n1 <- sum(val_label == classes[1])
  val_n2 <- sum(val_label == classes[2])
  val_n <- val_n1 + val_n2
  tropsvm.out <- tropsvm(data, label, assignment = best_assignment, ind = best_method_ind)
  cv.tropsvm.out <- list("coef" = tropsvm.out$coef,
                         "assignment" = best_assignment,
                         "method index" = best_method_ind,
                         "levels" = as.factor(classes),
                         "accuracy" = best_accuracy,
                         "nfold" = nfold)
  class(cv.tropsvm.out) <- "cv.tropsvm"
  cv.tropsvm.out
}

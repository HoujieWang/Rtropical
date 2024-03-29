% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.tropsvm.R
\name{predict.tropsvm}
\alias{predict.tropsvm}
\title{Predict Method for Tropical Support Vector Machines}
\usage{
\method{predict}{tropsvm}(object, newx, ...)
}
\arguments{
\item{object}{a fitted \code{tropsvm} object.}

\item{newx}{a data matrix, of dimension nobs x nvars used as testing data.}

\item{\dots}{Not used. Other arguments to predict.}
}
\value{
A vector of predicted values of a vector of labels.
}
\description{
Predicts values based upon a model trained by \code{tropsvm}.
}
\examples{

# data generation
library(Rfast)
e <- 100
n <- 10
N <- 10
s <- 5
x <- rbind(
  rmvnorm(n, mu = c(5, -5, rep(0, e - 2)), sigma = diag(s, e)),
  rmvnorm(n, mu = c(-5, 5, rep(0, e - 2)), sigma = diag(s, e))
)
y <- as.factor(c(rep(1, n), rep(2, n)))
newx <- rbind(
  rmvnorm(N, mu = c(5, -5, rep(0, e - 2)), sigma = diag(s, e)),
  rmvnorm(N, mu = c(-5, 5, rep(0, e - 2)), sigma = diag(s, e))
)
newy <- as.factor(rep(c(1, 2), each = N))

# train the tropical svm
tropsvm_fit <- tropsvm(x, y, auto.assignment = TRUE, ind = 1)

# test with new data
pred <- predict(tropsvm_fit, newx)

# check with accuracy
table(pred, newy)

# compute testing accuracy
sum(pred == newy) / length(newy)
}
\seealso{
\code{summary}, \code{coef} and the \code{tropsvm} function.
}

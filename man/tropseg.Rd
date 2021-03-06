% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tropseg.R
\name{tropseg}
\alias{tropseg}
\title{Compute Tropical Line Segment}
\usage{
tropseg(D1, D2, flag = FALSE)
}
\arguments{
\item{D1}{a vector of length 2.}

\item{D2}{a vector of length 2.}

\item{flag}{a logical value indicating if to swap \code{D1} and \code{D2}. (default: FALSE)}
}
\value{
A list containing the following:
\item{1}{End points as columns of each line component.}
\item{2}{Overall length of the tropical line segment}
}
\description{
Compute Tropical Line Segment
}
\references{
Lin, B., Monod, A. and Yoshida, R. (2020)
\emph{Tropical Geometric Variation of Phylogenetic Tree Shapes}
\url{https://arxiv.org/pdf/2010.06158.pdf}
}
\author{
Houjie Wang
Maintainer: Houjie Wang \email{wanghoujie6688@gmail.com}
}
\keyword{internal}

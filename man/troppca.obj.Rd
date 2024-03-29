% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/troppca.obj.R
\name{troppca.obj}
\alias{troppca.obj}
\alias{troppca.obj2}
\title{Compute Tropical PCA Objective}
\usage{
troppca.obj(pc, x_list, cl)

troppca.obj2(pc, x_list, cl)
}
\arguments{
\item{pc}{a matrix of principal components}

\item{x_list}{a list of vectors}

\item{cl}{cluster for parallel computing}
}
\value{
The numeric value of the objective function of tropical principle component analysis.
This is the sum of all tropical distance from each point to its projection on the tropical polytope.
}
\description{
Compute Tropical PCA Objective
}
\keyword{internal}

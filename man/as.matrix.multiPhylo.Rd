% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.matrix.multiPhylo.R
\name{as.matrix.multiPhylo}
\alias{as.matrix.multiPhylo}
\title{Vectorize a Set of Phylognetic Trees}
\usage{
\method{as.matrix}{multiPhylo}(x, tipOrder = x[[1]]$tip.label, parallel = FALSE, ncores = 2, ...)
}
\arguments{
\item{x}{an object of class \code{multiPhylo}}

\item{tipOrder}{a numeric vector of order of leaf names to which all trees
in the \code{multiPhylo} object will unified. If not specified on purpose,
the tip order of the first tree will be used.}

\item{parallel}{a logical value indicating if parallel computing should be used. (default: FALSE)}

\item{ncores}{a numeric value indicating the number of threads
utilized for multi-cored CPUs. (default: 2)}

\item{\dots}{Not used. Other arguments to as.vector}
}
\value{
A data matrix with each row a vector representation of a chronogram. Each element of the vector is the distance between two leaves.
}
\description{
Unifies tip labels of all phylogenetic trees in \code{multiPhylo}
object the same as the first tree and returns the
cophenetic distance of their corresponding chronogram.
}
\examples{
data(apicomplexa)
data <- as.matrix(apicomplexa[1: 10]) # matrixize first ten trees

}

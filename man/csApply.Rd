\name{csApply}
\alias{csApply}
\title{Apply a function over the intensities in a cytoSet}
\description{This a wrapper for \code{\link{sapply}} for
  objects of class \code{\link{cytoSet}}}
\usage{
csApply(X, FUN, ...)
}
\arguments{
  \item{X}{cytoSet.}
  \item{FUN}{the function to be applied.}
  \item{...}{optional arguments to \code{FUN}.}
}

\details{A wrapper for \code{\link{sapply}}.}

\value{Like \code{\link{sapply}}:
  If \code{FUN} always returns a scalar, then the value of
  this function is a named vector.
  If \code{FUN} always returns a vector of length \code{n},
  then the value of this function is an \code{n x length(X)}
  matrix with dimnames.
  Else, the value of this function is a named list whose
  value are the return values of the individual calls to \code{FUN}.}

\seealso{\code{\link{sapply}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/mga/whuber}}
\keyword{manip} 
\examples{
csApply(cset, nrow)
csApply(cset, colMeans)
}

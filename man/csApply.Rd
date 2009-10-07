\name{csApply}
\alias{csApply}
\title{Apply a function over the intensities in a cytoSet}
\description{This a wrapper for \code{\link[base:lapply]{sapply}} for
  objects of class \code{\link{cytoSet-class}}.}
\usage{
csApply(X, FUN, ..., simplify = TRUE)
}
\arguments{
  \item{X}{cytoSet.}
  \item{FUN}{the function to be applied.}
  \item{...}{optional arguments to \code{FUN}.}
  \item{simplify}{logical; should the result be simplified to a vector or
          matrix if possible? Gets passed on the \code{\link[base:lapply]{sapply}}.}
}

\details{A wrapper for \code{\link[base:lapply]{sapply}}.}

\value{Like \code{\link[base:lapply]{sapply}}:
  If \code{FUN} always returns a scalar, then the value of
  this function is a named vector.
  If \code{FUN} always returns a vector of length \code{n},
  then the value of this function is an \code{n x length(X)}
  matrix with dimnames.
  Else, the value of this function is a named list whose
  values are the return values of the individual calls to \code{FUN}.}

\seealso{\code{\link[base:lapply]{sapply}}}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\keyword{manip} 
\examples{
cset=readCytoSet(path=system.file("extdata", package="prada"),
  pattern="[A-Z][0-9][0-9]$")
csApply(cset, nrow)
csApply(cset, colMeans)
}

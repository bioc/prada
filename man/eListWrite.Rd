\name{eListWrite}
\alias{eListWrite}
\title{Write the eList}
\description{Writes the values of the eList of an eSet to a file}
\usage{
eListWrite(x,...)
}
\arguments{
  \item{x}{an \code{\link[Bioconductor:eSet-class]{eSet}}..}
  \item{...}{additional arguments passed to write.table.}
}

\details{The values of the eList of \code{x} are extracted and written
  to a file with write.table.}

\value{A an object of class}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\keyword{}
\examples{##see vignette 'How to apply the ddCt method'}

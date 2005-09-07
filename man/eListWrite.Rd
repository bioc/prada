\name{eListWrite}
\alias{eListWrite}
\title{Write the eList}
\description{Writes the values of the assayData list of an eSet to a file}
\usage{
eListWrite(x,...)
}
\arguments{
  \item{x}{an \code{\link[Bioconductor:eSet-class]{eSet}}..}
  \item{...}{additional arguments passed to write.table.}
}

\details{The values of the assayData list of \code{x} are extracted and written
  to a file with write.table. assayData of x must be a list.}

\value{A an object of class}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\keyword{}
\examples{##see vignette 'How to apply the ddCt method'}

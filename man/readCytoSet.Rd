\name{readCytoSet}
\alias{readCytoSet}
\title{Create a cytoSet object from one or more FCS 3.0 files}
\description{Create a cytoSet object from one or more FCS 3.0 files}
\usage{
  readCytoSet(files, path=".", ...)
}
\arguments{
  \item{files}{Optional character vector with filenames}
  \item{path}{Directory where to look for the files}
  \item{...}{If \code{files} is missing, the file names are obtained
     from \code{dir(path, ...)}}
}

\value{An object of class \code{link[cytoSet-class]{cytoSet}}.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\seealso{\code{\link{readFCSdata}}}

\examples{
readCytoSet(path=system.file("extdata", package="prada"), pattern="[A-Z][0-9][0-9]$")

}
\keyword{IO}


\name{readFCS}
\alias{readFCS}
\title{Read an FCS 3.0 file}
\description{Read an FCS 3.0 file: Data File Standard for Flow
  Cytometry, Version FCS3.0}
\usage{
  readFCS(filename)
}
\arguments{
  \item{filename}{Character of length 1: filename}
}

\details{For specifications of FACS 3.0 see
\url{http://www.isac-net.org} and the file
\url{../doc/fcs3.html} in the \code{doc} directory of the package.

The function works with the output of the FACS machine software from a
number of vendors. However, the FCS 3.0 standard includes some options
that are not yet implemented in this function. If you need extensions,
please let me know.
}

\value{An object of class \code{link[cytoFrame-class]{cytoFrame}}.}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{\code{\link{readFCSdata}}}

\examples{
sampdat <- readFCS(system.file("extdata", "fas Bcl2 plate323-04-04.A01",
                               package="prada"))
sampdat
exprs(sampdat[1:3,])
description(sampdat)[3:6]
}
\keyword{IO}


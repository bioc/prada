\name{readFCS}
\alias{readFCS}
\title{Read an FCS 3.0 file}
\description{Read an FCS 3.0 file: Data File Standard for Flow Cytometry, Version FCS3.0}
\usage{
  readFCS(filename)
}
\arguments{
  \item{filename}{Character of length 1: filename}
}

\details{For specifications of FACS 3.0 see
\url{http://www.isac-net.org} and the file
fcs3.html in the doc directory.

The FCS 3.0 has a many options that are not yet implemented in this
function. It works for the data produced by the FACS machine Calibur
from BD Bioscience together with the software Cellquest PRO, at the
Division of Molecular Genome Analysis at the German Cancer Research
Center. If you need extensions, please let me know.
}

\value{A matrix with column names, and atrribute "text".}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\seealso{}
\examples{
sampdat <- readFCS(system.file("extdata", "fas Bcl2 plate323-04-04.A01", package="prada"))
sampdat[1:5,]
dim(sampdat)
attr(sampdat, "text")[3:6]
}
\keyword{IO}


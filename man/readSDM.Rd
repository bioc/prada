\name{readSDM}
\alias{readSDM}
\title{Read an SDM file}
\description{Read an SDM file: Data Output File for SDS, Version 2.1}
\usage{
  readSDM(files, withoutPath=TRUE)
}
\arguments{
  \item{files}{Character vector of filenames}
  \item{withoutPath}{logical - if \code{TRUE}, the path to the files is
    stripped from the \code{filename} column of the resulting data frame.} 
}

\details{This function reads the data given in the colums
  'Detector','Sample' and 'Ct' of the specified SDM output files and
  stores them as a data.frame. An additional column including the
  respective filename is added.
}

\value{A data.frame}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}

\examples{
sampdat <- readSDM(system.file("extdata", "Experiment1.txt", package="prada"))
sampdat[1:5,]
}
\keyword{}


\name{statWell}
\alias{statWellLoess}
\alias{statWellParabola}
\title{Calculate a statistic for the data from one well}
\description{Calculate a statistic for the data from one well}
\usage{
statWellParabola(x, plotfile=as.character(NA), crosstalk)
statWellLoess(x, plotfile=as.character(NA), crosstalk, xmax, span) 
}
\arguments{
  \item{x}{data.frame with columns \code{brdu}, \code{trsf},
    \code{dapi}, \code{Field}, \code{cloneId}}
  \item{plotfile}{character. If not NA, a suitable plot is
    generated and saved as a PNG file under that name}
  \item{crosstalk}{numeric.
    Crosstalk factor between DAPI and CFP/YFP}
}

\details{}

\value{?}
\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}


\name{statWell}
\alias{statWell}
\title{Calculate a statistic for the data from one well}
\description{Calculate a statistic for the data from one well}
\usage{statWell(x, pars, plotfile=as.character(NA))}
\arguments{
  \item{x}{data.frame with columns \code{brdu}, \code{trsf},
    \code{dapi}, \code{Field}, \code{cloneId}}
  \item{pars}{list containing parameters for the calculation of the
    statistic. It must contain the following elements: \code{minRgTau},
    \code{minNrCells}, \code{ct}}
  \item{plotfile}{character. If not NA, a suitable plot is
    generated and saved as a PNG file under that name}
}

\details{}

\value{?}
\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}


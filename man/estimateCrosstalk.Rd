\name{estimateCrosstalk}
\alias{estimateCrosstalkWell}
\alias{estimateCrosstalkPlate}
\title{Estimate crosstalk between two color channels.}
\description{Estimate crosstalk between two color channels.}
\usage{
estimateCrosstalkPlate(x, plotfileprefix=NULL)
estimateCrosstalkWell(x, doPlot=TRUE)
}
\arguments{
  \item{x}{data frame.}
  \item{plotfileprefix}{character of length 1.}
  \item{doPlot}{logical of length 1.}
}

\details{This function is quite ideosyncratic - it was written for one
particular data sets and ought not be useful for any other purpose.}

\value{A numeric vector.}

\seealso{}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\keyword{manip}
\examples{##see vignette}

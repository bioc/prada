\name{devRes}
\alias{devRes}
\title{Resolution of current plotting device}
\description{Calculates what R thinks to be the resolution of the
  current graphic device.}
\usage{
devRes()
}

\details{
  This function may be used to get the resolution of the current
  graphics device. This can be important when calculating pixel
  coordinates for the output graphic. 
}

\value{A vector with items \code{xres} and
  \code{yres}, the resolution in \code{x} and \code{y} direction
respectively.}

\seealso{\code{\link{plotPlate}}}


\author{Florian Hahne}
\keyword{}
\examples{
devRes()
}

\name{vpLocation}
\alias{vpLocation}
\title{Absolute location of current viewport}
\description{Calculates the absolute location and size of the current grid
  viewport in inches an pixels.}

\usage{
vpLocation()
}

\details{
  This function may be used to get the absolte location of the current
  viewport on the current graphics device. It uses function
  \code{\link[prada::devRes]{devRes}} to get the device resolution for
  calculating pixel values. Locations are given by the two extreme
  coordinates in \code{x} and \code{y} direction.  
}

\value{A list with items \code{location, size, ilocation} and
  \code{isize}, the location and size of the viewport in pixels and
  iches respectively. }

\seealso{\code{\link[prada::plotPlate]{plotPlate}},
         \code{\link[prada::devRes]{devRes}}}


\author{Florian Hahne}
\keyword{}
\examples{
vpLocation()
}

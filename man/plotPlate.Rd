\name{plotPlate}
\alias{plotPlate}
\title{Plot a well statistic in false color representation.}
\description{Plot a well statistic in false color representation. The
  plot is supposed to resemble the physical geometry of a microtitre plate.}
\usage{
plotPlate(z, well, sizex=12, sizey=8, main, zlim=NULL, col, device,
          plotfile=NULL, width, na.action="zero")
}
\arguments{
  \item{z}{Numeric vector}
  \item{well}{Numeric vector of the same length as z. Its elements
    are expected to be integer in the range \code{1:(sizex*sizey)}.}
  \item{sizex}{Integer. The number of columns of the plate.}
  \item{sizey}{Integer. The number of rows of the plate.}
  \item{main}{Character of length 1. Plot title.}
  \item{zlim}{Range of z that is mapped into the color scale.}
  \item{col}{Character vector. Usually the names of 2 or 3 colors
    between which the color map is interpolated, using the function
    \code{\link{colorramp}}.}
  \item{device}{Character. Name of a valid graphics device,
    e.g. \code{x11}, \code{png}, \code{pdf}.}
  \item{plotfile}{The \code{file} argument of the call to the
    \code{device}, i.e. the name of the graphics output file.}
  \item{width}{The \code{width} argument of the call to the
    \code{device}, i.e. the width of the graphic.}
  \item{na.action}{Character. One of \code{"zero"} or \code{"omit"}.
    How should the wells for which \code{z} is NA be treated?
    For \code{"zero"}, they are plotted as if the value were 0. 
    For \code{"omit"}, they are omitted. 
}

\details{To do.}
\value{To do..}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{hplot}
\example{

}


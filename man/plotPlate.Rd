\name{plotPlate}
\alias{plotPlate}
\title{Plot a well statistic in false color representation.}
\description{Plot a well statistic in false color representation. The
  plot is supposed to resemble the physical geometry of a microtitre plate.}
\usage{
plotPlate(x, nrow=8, ncol=12, ind, main, xrange, col, 
  device, file, width, na.action="zero", desc)
}
\arguments{
  \item{x}{Numeric vector of length \code{ncol*nrow}}
  \item{nrow}{Numeric of length 1. The number of rows of the plate.}
  \item{ncol}{Numeric of length 1. The number of columns of the plate.}
  \item{ind}{Optional integer vector of equal length as \code{x}. It
    indicates the position of the respective value of
    x on the plate. Can be used to adress the problem of missing
    values. Each well that is not allocated a value of \code{x} by
    \code{ind} is set to \code{NA}.}
  \item{main}{Character of length 1. Plot title. This arguments gets passed 
     on to \code{\link[graphics:text]{text}} as the argument \code{labels}.}
  \item{xrange}{Range of \code{x} that is mapped into the color scale.}
  \item{col}{Character vector. Usually the names of two or three colors
    between which the color map is interpolated, using the function
    \code{\link[geneplotter::colorRampPalette]{colorRampPalette}}.}
  \item{device}{Character. Name of a valid graphics device,
    e.g. \code{x11}, \code{png}, \code{pdf}.}
  \item{file}{The \code{file} argument of the call to
    \code{device}, i.e. the name of the graphics output file.}
  \item{width}{The \code{width} argument of the call to 
    \code{device}, i.e. the width of the graphic.}
  \item{na.action}{Character. One of \code{"zero"} or \code{"omit"}.
    How should the wells for which \code{x} is NA be treated?
    For \code{"zero"}, they are plotted as if the value were 0. 
    For \code{"omit"}, they are omitted. }
  \item{desc}{Character of length 2. Legend for the two ends of the colorbar, e.g. 'act' and 'inh'. }
}

\details{}

\value{The function produces a plot in the graphics device specified by the
\code{device} argument.  

It returns a list with two elements. The element \code{which} is a vector with
the indices of those elements in \code{x} that were plotted (see argument
\code{na.action}). The element \code{coord} is a \code{length(which)} by 4
matrix in which each row specifies the corners of a rectangle that contains a
well. It is intended to be use as an argument to a subsequent call to
\code{\link{imageMap}}.}

\seealso{\code{\link[prada::imageMap]{imageMap}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{hplot}
\examples{
  plotPlate(runif(96), main="example 1", col=c("#0000e0", "#e00000"), width=7)
  plotPlate(runif(384), nrow=16, ncol=24,  
     main="example 2", col=c("#0000e0", "white", "#e00000"), width=7)
}


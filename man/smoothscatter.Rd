\name{smoothScatter}
\alias{smoothScatter}
\title{Smoothed scatterplots densities}
\description{Produces an image map showing empirical distribution of
  dots in scatterplots derived from two-dimensional histograms.
}
\usage{smoothScatter(x, y=NA, nrpoints=100, colrange=c("white",
"#231d1d"), nbin=200, bandwidth=20, ...) 
}
\arguments{
  \item{x}{Numeric vector containing x-values or n by 2 matrix containing x
    and y values.}
  \item{y}{Numeric vector containing y-values (optional). The length of \code{x}
    must be the same as that of \code{y}.}
  \item{nrpoints}{Numeric vector of length 1 giving number of points to
    be superimposed on image map.}
  \item{colrange}{Character vector with valid color
    identifiers (eg name or RGB values) from which a smooth color
    palette is derived.}
  \item{nbin}{Numeric vector containing the number of equally spaced points in each direction
    over which the density is to be estimated. Gets passed on to function
    \code{\link[KernSmooth::bkde2D]{bkde2D}}.}
  \item{bandwidth}{Numeric vector giving amount of
    smoothing in each coordinate direction.
    Gets passed on to function \code{\link[KernSmooth::bkde2D]{bkde2D}}.}
  \item{...}{Further arguments that are passed on to \code{image}.}
}

\details{Function is called for its effect of producing color images
  showing densities of scatterplots derived from smoothed two dimensional
  histograms (\code{\link[prada::smoothdens]{smoothdens}}. The range of the
  color palette used for the image can be specified and variable numbers of
  points can be superimposed on the image to show the spread of the
  data. If all points are to be plotted choose \code{nrpoints = Inf}. 
}

\value{
}  

\seealso{\code{\link[KernSmooth::bkde2D]{bkde2D}}}
\author{Florian Hahne}
\keyword{}
\examples{
sampdat <- readFCS(system.file("extdata", "fas Bcl2 plate323-04-04.A01",
package="prada"))
xxx <- exprs(sampdat[,1:2])
layout(matrix(1:4, ncol=2, byrow=TRUE))
smoothScatter(xxx, nrpoints=0)
smoothScatter(xxx)
smoothScatter(xxx, nrpoints=100, colrange=c("white", "blue"))
smoothScatter(xxx, nrpoints=Inf, colrange=c("white", "blue", "red"),
bandwidth=40)  
}


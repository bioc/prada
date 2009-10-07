\name{devDims}
\alias{devDims}
\title{Device Dimensions for plate plots}
\description{Calculate device dimensions for plate plots}
\usage{
devDims(width, height, ncol=12, nrow=8, res=72)
}
\arguments{
  \item{width}{Device width in inches.}
  \item{height}{Device width in inches.}
  \item{ncol}{Number of columns for plate plot.}
  \item{nrow}{Number of rows for plate plot.}
  \item{res}{The resolution of the graphic device used for plotting.}
}

\details{
  The function computes the device dimensions needed to create plate
  plots that fit perfectly in the device. This is necessary to retain
  the aspect ratio of the plots.

  One of width or height need to be specified, the missing value will be
  computed.
}

\value{A list with items \code{width, height, pwidth} and
  \code{pheight}. These are the width and height values in inches and
  pixels respectively.}

\seealso{\code{\link{plotPlate}}}


\author{Florian Hahne}
\keyword{}
\examples{
devDims(width=10)
}

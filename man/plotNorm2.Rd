\name{plotNorm2}
\alias{plotNorm2}
\title{Plot fitted bivariate normal distribution.}
\description{Plots objects derived from function
  \code{\link[prada::fitNorm2]{fitNorm2}} in false color represenation.
}
\usage{
plotNorm2(fn, colrange=c("gray82", "blue"), center=TRUE, selection=FALSE,
                      ellipse=FALSE, pch=20, cex=1, col="dens", ...)
}
\arguments{
  \item{fn}{List. Object derived from  function
    \code{\link[prada::fitNorm2]{fitNorm2}}}
  \item{colrange}{Character vector with valid color
    identifiers (eg name or RGB values) from which a smooth color
    palette is derived.}
  \item{center}{Logical. Assign center of distribution.}
  \item{selection}{Logical. Mark all points beyond selection.}
  \item{ellipse}{Logical. Plot area and borders of selection as ellipse.}
  \item{pch}{see \code{par}}
  \item{cex}{see \code{par}}
  \item{col}{see \code{par} or special cases \code{dens} for coloring
    according to density and \code{prob} for coloring according to
    probability.}
  \item{...}{further arguments that are passed on to \code{plot}.}
}

\details{Produces a scatterplot of paired data showing the densities of
  the fitted bivariate distribution from function
  \code{\link[prada::fitNorm2]{fitNorm}} in false color
  representation. Additionally a selection of data points can be
  highlighted either by marking outliers or by showing its area.
}

\value{A list containing items \code{p, cov, mu, S} (density values for
  each data pair, resulting object from call to cov.rob, midpoint of
  distribution, covariance matrix).
}  

\seealso{\code{\link[prada::fitNorm2]{fitNorm2}}, \code{\link[prada::addEllipse]{addEllipse}}}

\author{Florian Hahne}
\keyword{}
\examples{
 sampdat <- readFCS(system.file("extdata",
   "fas-Bcl2-plate323-04-04.A01", package="prada"))
 nfit   <- fitNorm2(exprs(sampdat[,1:2]), scalefac=2)
 plotNorm2(nfit, selection=TRUE, ellipse=TRUE)
}
 


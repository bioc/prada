\name{fitNorm2}
\alias{fitNorm2}
\title{Fit bivariate normal distribution.}
\description{Fits a bivariate normal distribution into a data set of
  paired values and selects data points according to their standard
  deviation from the fitted distribution.
}
\usage{fitNorm2(x, y=NA, scalefac=1, method="covMcd")
}
\arguments{
  \item{x}{Numeric vector containing x-value or n by 2 matrix containing x and
    y values.}
  \item{y}{Numeric vector containing y-value (optional). The length of \code{x}
    must be the same as that of \code{y}.}
  \item{scalefac}{Numeric vector giving factor of standard deviations
    used for data selection (all points within \code{scalefac} standard
    deviations are selected).}
  \item{method}{One of \code{covMcd} or \code{cov.rob} defining method
    used for computation of covariance matrix.}
}

\details{Computes the densities of a bivariate normal distribution from
  the covariance matrix of the paired data. Covariance matrices are
  acquired either by function \code{\link[rrcov::covMcd]{covMcd}}
  (considerably faster) or by function \code{\link[MASS::cov.rob]{cov.rob}}.
}

\value{A list containing items \code{mu} (midpoint of distribution),
  \code{S} (covariance matrix), \code{p} (density values for each
  data pair), \code{sel} (selection of data points), \code{scalefac}
  (factor of standard deviations used for data selection), \code{data}
  (x and y values of data points).}
}  

\seealso{\code{\link[MASS::cov.rob]{cov.rob}},
  \code{\link[rrcov::covMcd]{covMcd}},
  \code{\link[prada::plotNorm2]{plotNorm2}}}

\author{Florian Hahne}
\keyword{}
\examples{
 sampdat <- readFCS(system.file("extdata",
   "fas Bcl2 plate323-04-04.A01", package="prada"))
 nfit    <- fitNorm2(exprs(sampdat[,1:2]), scalefac=2)
 plotNorm2(nfit, selection=TRUE, ellipse=TRUE) 
}


\name{threshold}
\alias{threshold}
\title{Discretize a two-dimensional data space into quadrants by
  applying thresholds}

\description{Discretize a two-dimensional data space into quadrants by
  applying thresholds.}

\usage{threshold(x, y, xthr, ythr)}

\arguments{
  \item{x}{Vector containing x or matrix containing x and y values of
     bivariate data.}
  \item{y}{Optional vector containing y values of bivariate data.}
  \item{xthr}{x value seperating 'left' and 'right'.}
  \item{ythr}{y value seperating `up` and 'down'.}
}

\details{The function returns a 2x2 matrix giving the counts
for each quadrant. Events with values equal to the thresholds are
counted to the left or down respectively.}

\value{2x2 matrix.}  

\seealso{}

\author{Florian Hahne}
\keyword{}
\examples{
threshold(cbind(c(1, 1, 2, 2, 2, 4), c(1, 4, 2, 4, 5, 4)), xthr=3, ythr=3)
}


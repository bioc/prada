\name{colorramp}
\alias{colorramp}
\title{Interpolate between colors}
\description{Interpolate between colors}
\usage{colorramp(col)}
\arguments{
  \item{col}{A vector of colors, as in \code{\link[graphics:col2rgb]{col2rgb}}.}
}

\details{The function is adapated from the function
  \code{\link[geneplotter:colorRamp]{colorRamp}} in the package
  \code{geneplotter}.}

\value{A function that takes an integer and produces
  a color map which interpolates between the colors in \code{col}.}

\author{Thomas Lumley}

\examples{
  cr <- colorramp(c("green", "red"))
  cr(5)
}
\keyword{arith}


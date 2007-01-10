\name{barploterrbar}
\alias{barploterrbar}
\title{Barplot with error bars.}
\description{Barplot with error bars.}
\usage{
barploterrbar(y, yl, yh, barcol="orange", errcol="black", horiz=FALSE,
w=0.2, ylim=c(0, max(yh)*1.05), ...)
}
\arguments{
  \item{y}{Numeric vector.}
  \item{yl}{Numeric vector of same length as y.}
  \item{yh}{Numeric vector of same length as y.}
  \item{barcol}{Color of the bars.}
  \item{errcol}{Color of the error bars.}
  \item{horiz}{Logical. As in \code{\link[graphics:barplot]{barplot}}.}
  \item{w}{The plot limits. The default value wil cause the error bars
    to fit nicely on the plotting device.}
  \item{ylim}{Size of the error bar ticks.}
  \item{...}{Further arguments that get passed on to
     \code{\link[graphics:barplot]{barplot}}.}
}

\details{The function calls \code{\link[graphics:barplot]{barplot}} with
\code{y} and decorates it with error bars according to \code{yl} and
\code{yh}.}

\value{The function is called for its side effect, producing a plot.} 

\seealso{\code{\link[graphics:barplot]{barplot}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{hplot}
\examples{
  y <- matrix(runif(80), ncol=5)
  ym <- apply(y, 2, mean)
  dy <- apply(y, 2, sd)*2/sqrt(nrow(y))
  barploterrbar(ym, ym-dy, ym+dy, barcol="#0000c0", errcol="orange",
    ylim=c(0, max(ym+dy)))
}


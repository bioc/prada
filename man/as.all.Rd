\name{as.all}
\alias{as.all}
\title{Coercion without introduction of NAs}
\description{Coercion without introduction of NAs}
\usage{
as.all(x, what)
}
\arguments{
  \item{x}{an object.}
  \item{what}{character of length 1.} 
}

\details{
  The function calls
  \code{do.call(paste("as.", what, sep=""), list(x))},
  and checks whether any \code{NAs} were introduced.
}

\value{A vector of type \code{what}}

\seealso{\code{\link[methods:as]{as}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}
\examples{
  as.all(runif(5)*10, "integer")
}

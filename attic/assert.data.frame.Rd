\name{assert.data.frame}
\alias{assert.data.frame}
\title{Check whether the argument is a data frame}
\description{Check whether the argument is a data frame}

\usage{
  assert.data.frame(x, onerow=TRUE)
}

\arguments{
  \item{x}{An S object.}
  \item{onerow}{Logical: check whether x (if it is a data frame) has
    exactly one row.}
}

\details{
  This function is intended to be called subsequently to calls to
  \code{\link{sqlQuery}} (from package RODB) that contain a SELECT
  statement.
  If the query could be answered successfully,
  \code{\link{sqlQuery}} returns with a data.frame. Else, it returns 
  a character vector with the error message. Thus,
  \code{assert.data.frame} checks whether \code{x} is a data frame,
  returns TRUE if yes, and calls \code{\link{stop}} if no.
}

\value{  
  If x is a data.frame (with one row), the function returns TRUE.
}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}

\seealso{\code{\link{sqlQuery}},\code{\link{stop}},\code{\link{stopifnot}}}

\examples{
  x = data.frame(runif(100), runif(100))
  assert.data.frame(x, onerow=FALSE)
}

\keyword{programming}
\keyword{error}

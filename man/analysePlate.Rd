\name{analysePlate}
\alias{analysePlate}
\title{Apply a statistic to the data from each well in a plate}
\description{Apply a statistic to the data from each well in a plate}
\usage{
analysePlate(x, wellcol="well", statfun, platename, plotdir=".", ...)
}
\arguments{
  \item{x}{data frame. It must contain a column whose name is the value
    of \code{wellcol}, and further columns that are needed by the function
    named by \code{stat}.}
  \item{wellcol}{character of length 1. Name of a column in \code{x} that
    contains the well ID. The data type of the column can be character,
    integer, or factor.}
  \item{statfun}{character of length 1. Name of a function that can calculate
    a statistic from selected rows of \code{x}.}
  \item{platename}{character of length 1. The name or ID of this plate, which
  will be used for graphics output filenames and as the value of the column
  \code{platename} of the return value.}
  \item{plotdir}{character of length 1. The name of directory where
    diagnostic plots will be saved.}
  \item{...}{further arguments that are passed on to \code{statfun}.}
}

\details{The semantics of this function are similar to
\code{\link[base:tapply]{tapply}}, but some additional checking and 
reporting is performed, and the return value is a data frame.}

\value{A data frame with as many rows as there are unique elements of 
\code{x[, wellcol]}. It has a (among others) a column named \code{platename}
whose value in each row is the value of the corresponding input argument.}

\seealso{\code{\link{statWellLocfit}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}
\examples{##see vignette}

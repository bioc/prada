\name{readFCSaux}
\alias{readFCSdata}
\alias{readFCStext}
\alias{readFCSheader}
\alias{readFCSgetPar}
\title{Auxiliary functions for readFCS}
\description{Auxiliary functions for readFCS - not normally called by
  the user}
\usage{
  readFCSgetPar(x, pnam)
  readFCSheader(con)
  readFCStext(con, offsets)
  readFCSdata(con, offsets, x)
}
\arguments{
  \item{x}{Named character vector.}
  \item{pnam}{Character vector, its elements must be contained in
    \code{names(x)}.}
  \item{con}{Connection.}
  \item{offsets}{Integer vector of length 6 with byte offsets of the
    header, text, and data blocks.}
}

\details{These functions are not normally called by the user. See
  \code{\link{readFCS}} instead.}

\value{Various.}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{\code{\link{readFCS}}}
\examples{}
\keyword{IO}


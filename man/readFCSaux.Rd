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
  readFCStext(con, offsets, delimiter="\\")
  readFCSdata(con, offsets, x, endian="big")
}
\arguments{
  \item{x}{Named character vector.}
  \item{pnam}{Character vector, its elements must be contained in
    \code{names(x)}.}
  \item{con}{Connection.}
  \item{offsets}{Integer vector of length 6 with byte offsets of the
    header, text, and data blocks.}
  \item{delimiter}{Character.}
  \item{endian}{Character. Specifies the byte order, see \code{\link[base:readBin]{readBin}}.}
}

\details{These functions are not normally called by the user.}

\value{Various.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\seealso{\code{\link{readFCS}}}
\examples{}
\keyword{IO}


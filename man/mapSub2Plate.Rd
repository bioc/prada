\name{mapSub2Plate}
\alias{mapSub2Plate}
\title{Map subplate well ID and subplate ID to 96-well plate well ID}
\description{Map subplate well ID (e.g. from chamberslides) and subplate ID to 96-plate well ID}
\usage{mapSub2Plate(sub, well)}
\arguments{
  \item{sub}{Character. Subplate identifier, see details.}
  \item{well}{Integer. Well identifier, see details.}
}
\details{\code{sub} and \code{well} must have the same length.

  For chamberslides, \code{sub} has the values
  "a", "b", "c12", or "c34". Alternatively, \code{sub} can have the values
  "1-48", "49-96".

  \code{well} must be an integer number in the appropriate range.
}

\value{An numeric vector of the same length as \code{well} with the
  plate well IDs.
}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\seealso{}
\examples{
  plcat = function(iw, x) {
    d = numeric(96)
    d[iw] = x
    for (i in 0:7)
      cat(sapply(d[(1:12)+i*12], function(s) { sprintf('\%4s', s) }), '\n')
    cat('\n')
  }

  sub  = c(rep("a", 32), rep("b", 32), rep("c12", 16), rep("c34", 16))
  well = rep(1:32, 3)
  iw   = mapSub2Plate(sub, well)
  plcat(iw, paste(well))
  plcat(iw, sub)

  sub  = c(rep("1-48",  48), rep("49-96", 48))
  well = 1:96
  iw   = mapSub2Plate(sub, well)
  plcat(iw, paste(well))
  plcat(iw, sub)
}

\keyword{manip}


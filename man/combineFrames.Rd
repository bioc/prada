\name{combineFrames}
\alias{combineFrames}

\title{Combine the cytoFrames within a cytoSet according to some
  grouping factor}
\description{Combine the cytoFrames within a cytoSet according to some
  grouping factor.}

\usage{combineFrames(x, by)}
\arguments{
\item{x}{cytoSet.}
\item{by}{factor. Length must be same as that of \code{x}.}
}


\value{cytoSet.}  


\author{Wolfgang Huber <huber@ebi.ac.uk>}

\examples{
cset <- readCytoSet(path=system.file("extdata", package="prada"),
                 pattern="[A-Z][0-9][0-9]$")
nr1 <- csApply(cset, nrow)
sm1 <- csApply(cset, sum)

fac <- factor(c(1,1,2,2,2,2))
cc  <- combineFrames(cset, fac)
nr2 <- csApply(cc, nrow)
sm2 <- csApply(cc, sum)

stopifnot(all(nr2==tapply(nr1, fac, sum)))
stopifnot(all(sm2==tapply(sm1, fac, sum)))

}


\name{getDye}
\alias{getDye}
\title{Extract the Dye name (CFP or YFP) from the clone ID by pattern matching}
\description{Extract the Dye name (CFP or YFP) from the clone ID by pattern matching}
\usage{getDye(x)}
\arguments{\item{x}{Character vector with clone ID(s)}}

\details{The following patterns are recognised:
  \enumerate{
    \item ...[CY], ...[CY] Numbers, ...[CY] Numbers Letters, 
    \item [CY]...
    \item letter letter [CY]...
  }
See also the source code.}

\value{Charatcer vector of the same length as \code{x}.}
\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\examples{getDye("dorit10_3c_1-2--01-pDhfbr2_78k24_1E3Y1")}
\keyword{manip}


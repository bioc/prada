\name{removeCensored}
\alias{removeCensored}
\alias{removeCensored,matrix-method}
\alias{removeCensored,data.frame-method}

\title{Remove rows that contain censored data}
\description{Remove rows that contain censored data in the
  colums of \code{x} specified by \code{columns}.}

\usage{
\S4method{removeCensored}{matrix}(x, values, columns, na.rm=TRUE)
\S4method{removeCensored}{data.frame}(x, values, columns, na.rm=TRUE)
\S4method{removeCensored}{cytoFrame}(x, values, columns, na.rm=TRUE)
}
\arguments{
\item{x}{Object of class matrix, data.frame, or cytoFrame.}
\item{values}{Values that correspond to censored data. If missing,
  \code{range(x)} is used.}
\item{columns}{Numeric or character vector specifying the columns
  of \code{x} that are compared against \code{values}. If missing,
  \code{1:ncol(x)} is used.}
\item{na.rm}{Logical. If \code{TRUE}, rows that contain \code{NA}
  values are also removed.}
}

\details{The function removes all rows that contain, in the columns
  specified by the \code{columns} argument, values that are contained in
  the \code{values} argument. If \code{na.rm} is \code{TRUE}, then
  rows that contain \code{NA} values are also removed.

  An application is with FACS data, where measurements outside of
  the detector's dynamic range produce minimal
  or maximal values. For example, if a 16-bit A/D converter was used,
  top-censored data would have a value of 65535.}

\value{Object of the same class as \code{x}, with some rows removed.}  

\seealso{}

\author{Florian Hahne <f.hahne@dkfz.de>, Wolfgang Huber <huber@ebi.ac.uk>}
\keyword{}
\examples{
set.seed(8215)
mat <- matrix(floor(runif(20000)*1024), ncol=4)
range(mat[,1])
mat <- removeCensored(mat, columns=1:2)
range(mat[,1])
range(mat[,3])
}


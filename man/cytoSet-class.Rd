\name{cytoSet-class}
\docType{class}

\alias{cytoSet-class}

\alias{[,cytoSet-method}
\alias{[[,cytoSet-method}

\alias{show,cytoSet-method}
\alias{length,cytoSet-method}
\alias{pData,cytoSet-method}

\alias{colnames,cytoSet-method}
\alias{colnames<-,cytoSet-method}

\alias{phenoData,cytoSet-method}
\alias{phenoData<-,cytoSet,phenoData-method}

\title{'cytoSet': a class for storing raw data from
a quantitative cell-based assay}

\description{This class is a container for a set of
  \code{\link[cytoFrame-class]{cytoFrame}} objects}

\usage{
\S4method{phenoData}{cytoSet}(object)
\S4method{phenoData<-}{cytoSet,ANY}(object,value)
\S4method{pData}{cytoSet}(object)
\S4method{colnames}{cytoSet}(object)
\S4method{colnames<-}{cytoSet,ANY}(object,value)
\S4method{show}{cytoSet}(object)
\S4method{[}{cytoSet,ANY}(x, i)
\S4method{[[}{cytoSet,ANY}(x, i)
}

\arguments{
\item{object,x}{Objects of class \code{cytoSet}.}
\item{value}{Replacement value.}
\item{i}{Index.}
}

\section{Creating Objects}{
Objects can be created using the function
\code{\link{readCytoSet}} or via\cr
\code{  new('cytoSet,}\cr
\code{    frames = ....,   # environment with cytoFrames}\cr
\code{    phenoData = .... # object of class phenoData}\cr
\code{    colnames = ....  # object of class character}\cr
\code{  )}\cr
}

\section{Slots}{\describe{
\item{\code{frames}:}{An \code{\link[base:environment]{environment}}
containing one or more \code{\link[cytoFrame-class]{cytoFrame}} objects.}

\item{\code{phenoData}:}{A
\code{\link[Biobase:phenoData-class]{phenoData}}. Each row
corresponds to one of the cytoFrames in the \code{frames} slot.
It is mandatory that the pData has column named \code{framename}}

\item{\code{colnames}:}{A \code{character} object with the
(common) column names of all the data matrices in the cytoFrames.}
}}

\section{Methods}{\describe{
\item{[, [[}{subsetting. If \code{x} is \code{cytoSet}, then \code{x[i]}
returns a \code{cytoSet} object, and \code{x[[i]]} a \code{cytoFrame}
object. The semantics is similar to the behavior of the subsetting
operators for lists.}
\item{colnames, colnames<-}{extract or replace the \code{colnames} slot.}
\item{phenoData, phenoData<-}{extract or replace the \code{phenoData} slot.}
\item{show}{display summary.}
}}

\author{Wolfgang Huber}

\seealso{\code{\link{readCytoSet}}, \code{\link{cytoFrame-class}}}
\examples{
cset=readCytoSet(path=system.file("extdata", package="prada"),
  pattern="[A-Z][0-9][0-9]$")
cset
pData(cset)
cset[[1]]
cset[["fas Bcl2 plate323-04-04.A02"]]
cset["fas Bcl2 plate323-04-04.A02"]
cset[1:3]
}
\keyword{classes}


\name{cytoSet-class}
\docType{class}
\alias{cytoSet}
\alias{cytoSet-class}
\alias{[,cytoSet-method}
\alias{[[,cytoSet-method}
\alias{[[<-,cytoSet-method}
\alias{show,cytoSet-method}
\alias{length,cytoSet-method}
\alias{pData,cytoSet-method}
\alias{colnames,cytoSet-method}
\alias{colnames<-,cytoSet-method}
\alias{phenoData,cytoSet-method}
\alias{phenoData<-,cytoSet,phenoData-method}
\alias{plot,cytoSet,missing-method}


\title{'cytoSet': a class for storing raw data from
a quantitative cell-based assay}

\description{This class is a container for a set of
  \code{\link[prada]{cytoFrame}} objects}

%\usage{
%\S4method{phenoData}{cytoSet}(object)
%\S4method{phenoData<-}{cytoSet,ANY}(object,value)
%\S4method{pData}{cytoSet}(object)
%\S4method{colnames}{cytoSet}(object)
%\S4method{colnames<-}{cytoSet,ANY}(object,value)
%\S4method{show}{cytoSet}(object)
%\S4method{[}{cytoSet,ANY}(x, i)
%\S4method{[[}{cytoSet,ANY}(x, i)
%}

%\arguments{
%\item{object,x}{Objects of class \code{cytoSet}.}
%\item{value}{Replacement value.}
%\item{i}{Index.}
%}

\section{Creating Objects}{
Objects can be created using the function
\code{\link[prada]{readCytoSet}} or via\cr
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
It is mandatory that the pData has column named \code{name}}
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
\item{plot}{Scatterplot of one or all (consecutively) \code{\link[cytoFrame]{cytoFrame}}
objects.  The additional argument \code{gate} can be used to plot
subsets of the data defined by an object of class
\code{\link[prada]{gate}} or \code{\link[prada]{gateSet}}.}
}}

\section{Important note on storage and performance}{
The bulk of the data in a \code{cytoSet} object is stored in an
\code{\link[base:environment]{environment}}, and is therefore not
automatically copied when the \code{cytoSet} object is copied. If
\code{x} is an object of class \code{cytoSet}, then the code
\preformatted{y <- x} will create a an object \code{y} that contains
copies of the \code{\link[Biobase:phenoData-class]{phenoData}} and
administrative data in \code{x}, but refers to the \emph{same}
environment with the actual fluorescence data. See below for how to
create proper copies.

The reason for this is performance. The pass-by-value semantics of
function calls in \code{R} can result in numerous copies of the same
data object being made in the course of a series of nested function
calls. If the data object is large, this can result in a considerable
cost of memory and performance. \code{cytoSet} objects are intended to
contain experimental data in the order of hundreds of Megabytes,
which can effectively be treated as read-only: typical tasks
are the extraction of subsets and the calculation of summary statistics.
This is afforded by the design of the \code{cytoSet}
class: an object of that class contains a \code{phenoData} slot,
some administrative information, and a \emph{reference} to an
environment with the fluorescence data; when it is copied, only the
reference is copied, but not the potentially large set of fluorescence
data themselves.

However, note that subsetting operations, such as
\preformatted{y <- x[i]}
do create proper copies, including a copy of the appropriate part of the
fluorescence data, as it should be expected. Thus, to make a proper copy
of a \code{cytoSet} \code{x}, use
\preformatted{y <- x[seq(along=x)]}
}

\author{Florian Hahne, Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\seealso{\code{\link[prada]{readCytoSet}},
\code{\link[prada]{cytoFrame}}, \code{\link[prada]{gate}},
\code{\link[prada]{gateSet}}}

\examples{
cset<-readCytoSet(path=system.file("extdata", package="prada"),
  pattern="[A-Z][0-9][0-9]$")
cset
pData(cset)
cset[[1]]
cset[["fas-Bcl2-plate323-04-04.A02"]]
cset["fas-Bcl2-plate323-04-04.A02"]
cset[1:3]

cset[[1]] <- exprs(cset[[1]])[1:100, ]

\dontrun{
plot(cset[2])
}

}
\keyword{classes}


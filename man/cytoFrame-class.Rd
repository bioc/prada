\name{cytoFrame-class}
\docType{class}

\alias{cytoFrame-class}

\alias{[,cytoFrame-method}
% \alias{[<-,cytoFrame-method}

\alias{description,cytoFrame-method}
\alias{description<-,cytoFrame,character-method}

\alias{exprs,cytoFrame-method}
\alias{exprs<-,cytoFrame,matrix-method}

\alias{colnames,cytoFrame-method}
\alias{colnames<-,cytoFrame-method}

\alias{show,cytoFrame-method}

\alias{$.cytoFrame}

\title{'cytoFrame': a class for storing observed quantitative
  properties from a population of cells, for example from a FACS run or
  from automated microscopy}

\description{This class represents the data contained in a
  FCS 3.0 file or similar data structures.}

\usage{
\S4method{exprs}{cytoFrame}(object)
\S4method{exprs<-}{cytoFrame,matrix}(object,value)
\S4method{description}{cytoFrame}(object)
\S4method{description<-}{cytoFrame,character}(object,value)
\S4method{colnames}{cytoFrame}(object)
\S4method{colnames<-}{cytoFrame,ANY}(object,value)
\S4method{show}{cytoFrame}(object)
\S4method{[}{cytoFrame,ANY}(x, i, j, ..., drop=FALSE)
% \S4method{[<-}{cytoFrame,ANY,ANY,...,matrix}(x, i, j, ..., value)
$.cytoFrame(x, val)
}

\arguments{
\item{object,x}{Objects of class \code{cytoFrame}.}
\item{value}{Replacement value.}
\item{i,j}{Indices.}
\item{...}{Further arguments that get passed on.}
}


\section{Creating Objects}{
Objects can be created using\cr
\code{  new('cytoFrame,}\cr
\code{    exprs       = ...., # Object of class matrix}\cr
\code{    description = ....  # Object of class character}\cr
\code{  )}\cr
or the function \code{\link{readFCS}}.}

\section{Slots}{\describe{
\item{\code{exprs}:}{Object of class \code{matrix} containing
      the measured intensities. Rows correspond to cells, columns to the
      different channels. The \code{colnames} attribute of the matrix is
      supposed to hold the names or identifiers for the channels. The
      \code{rownames} attribute would usually not be set.}
\item{\code{description}:}{A named character vector
      containing the experiment description as key-value pairs.}
}}

\section{Methods}{\describe{
\item{[}{subsetting. Returns an object of class \code{cytoFrame}.
  The subsetting is applied to the \code{exprs} slot, while the
  \code{description} slot is unchanged.}
% \item{[<-}{replace subsets.}
\item{exprs, exprs<-}{extract or replace the intensities.}
\item{description, description<-}{extract or replace the description.}
\item{show}{display summary.}
}}

\details{
FCS 3.0 is the Data File Standard for Flow Cytometry, Version FCS3.0.
}
			
\author{Wolfgang Huber}

\seealso{\code{\link{readFCS}}, \code{\link{cytoSet-class}}}
\examples{

intens <- matrix(runif(100), ncol=4)
colnames(intens) <- c("FL1-H", "FL2-H", "FL3-H", "FL4-H")

a <- new("cytoFrame",
          exprs=intens,
          description=c(name="example data", date=date()))

description(a)
dim(exprs(a))           

a[1:3, -4]

}
\keyword{classes}


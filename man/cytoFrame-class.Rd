\name{cytoFrame-class}
\docType{class}
\alias{cytoFrame-class}
\alias{cytoFrame}
\alias{[,cytoFrame,ANY,ANY,ANY-method}
\alias{description,cytoFrame-method}
\alias{description<-,cytoFrame,character-method}
\alias{exprs,cytoFrame-method}
\alias{exprs<-,cytoFrame,matrix-method}
\alias{colnames,cytoFrame-method}
\alias{colnames<-,cytoFrame-method}
\alias{show,cytoFrame-method}
\alias{plot,cytoFrame,missing-method}
\alias{$.cytoFrame}
\alias{gate,cytoFrame-method}
\alias{gate<-,cytoFrame,gate-method}
\alias{gate<-,cytoFrame,gateSet-method}
\alias{drawGate,cytoFrame-method}
\alias{drawGate,cytoFrame-method}
\alias{drawGate,matrix-method}
\alias{appendGates,cytoFrame-method}
\alias{ncol,cytoFrame-method}
\alias{nrow,cytoFrame-method}
\alias{applyGate,cytoFrame,gate-method}
\alias{applyGate,cytoFrame,gateSet-method}
\alias{applyGate,cytoFrame,character-method}
\alias{applyGate,cytoFrame,numeric-method}
\alias{applyGate,cytoFrame,logical-method}
\alias{hist,cytoFrame-method}


\title{'cytoFrame': a class for storing observed quantitative
  properties from a population of cells, most likely from a FACS run or, 
  alternatively, from automated microscopy}

\description{This class represents the data contained in a
  FCS 3.0 file or similar data structures.}

%\usage{
%\S4method{exprs}{cytoFrame}(object)
%\S4method{exprs<-}{cytoFrame,matrix}(object,value)
%\S4method{description}{cytoFrame}(object)
%\S4method{description<-}{cytoFrame,character}(object,value)
%\S4method{colnames}{cytoFrame}(object)
%\S4method{colnames<-}{cytoFrame,ANY}(object,value)
%\S4method{show}{cytoFrame}(object)
%\S4method{[}{cytoFrame,ANY}(x, i, j, ..., drop=FALSE)
%\S4method{plot}{cytoFrame,missing}(x, ...)
% \S4method{[<-}{cytoFrame,ANY,ANY,...,matrix}(x, i, j, ..., value)
%$.cytoFrame(x, val)
%}

%\arguments{
%\item{object,x}{Objects of class \code{cytoFrame}.}
%\item{value}{Replacement value.}
%\item{i,j}{Indices.}
%\item{...}{Further arguments that get passed on.}
%}

\section{Creating Objects}{
Objects can be created using\cr
\code{  new('cytoFrame',}\cr
\code{    exprs       = ...., # Object of class matrix}\cr
\code{    description = ....  # Object of class character}\cr
\code{  )}\cr
or the function \code{\link[prada]{readFCS}}.}

\section{Slots}{\describe{
\item{\code{exprs}:}{Object of class \code{matrix} containing
      the measured intensities. Rows correspond to cells, columns to the
      different channels. The \code{colnames} attribute of the matrix is
      supposed to hold the names or identifiers for the channels. The
      \code{rownames} attribute would usually not be set.}
\item{\code{description}:}{A named character vector
      containing the experiment description as key-value pairs.}
\item{\code{well}:}{A single integer vector giving the position of the
     well on a microtitre plate. This only applies when using the object
     within a \code{\link[prada]{cytoSet}} collection and will usually
     be filled in by the function \code{\link[prada]{readCytoSet}}.}
\item{\code{gate}:}{An object of class
      \code{\link[prada]{gateSet}}. This object can be used to
      select defined subsets of the data, a process referred to as
      \code{gating} in the analysis of flow-cytometry data.}
}}

\section{Methods}{\describe{
\item{[}{subsetting. Returns an object of class \code{cytoFrame}.
  The subsetting is applied to the \code{exprs} slot, while the
  \code{description} slot is unchanged.}
\item{exprs, exprs<-}{extract or replace the intensities.}
\item{description, description<-}{extract or replace the description.}
\item{show}{display summary.}
\item{plot}{scatterplot for \code{cytoFrame} objects. The additional
argument \code{gate} can be used to plot subsets of the data defined by
either an object of class \code{\link[prada]{gate}} or by a character
vector giving the name of one of the gates in the list.}
\item{gate,gate<-}{extract or replace the list of gates.}
\item{ncol,nrow}{extract the dimensions of the data matrix.}
\item{appendGate}{Append a gate or gateSet to the gate slot.}
\item{drawGate}{Create an object of class \code{\link[prada]{gate}} or
\code{\link[prada]{gateSet}} based on a selection made from the
data. See \code{\link[prada]{gateMatrix}} for further details.}
\item{hist}{Draw a histogram of the data}
}}

\details{
Although objects of class \code{cytoFrame} can be used to hold arbitrary
data of cell populations, the main focus lies on flow-cytometry data.
 
FCS 3.0 is the Data File Standard for Flow Cytometry, Version
FCS3.0. See the vignette of this package for additional information on
using the object system for handling of flow-cytometry data. 
}
			
\author{Florian Hahne, Wolfgang Huber}

\seealso{\code{\link[prada]{readFCS}},
\code{\link[prada]{cytoSet}}, \code{\link[prada]{gate}},
\code{\link[prada]{gateSet}}, \code{\link[prada]{gateMatrix}}}

\examples{
intens <- matrix(runif(100), ncol=4)
colnames(intens) <- c("FL1-H", "FL2-H", "FL3-H", "FL4-H")

a <- new("cytoFrame",
          exprs=intens,
          description=c(name="example data", date=date()))

description(a)
dim(exprs(a))           

a[1:3, -4]

plot(a)
\dontrun{
g1 <- drawGate(a, name="Gate1")
}
}

\keyword{classes}


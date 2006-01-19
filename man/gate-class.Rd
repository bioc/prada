\name{gate-class}
\docType{class}
\alias{gate-class}
\alias{gate}
\alias{applyGate,gate,matrix-method}
\alias{applyGate,gate,cytoFrame-method}
\alias{show,gate-method}
\alias{names,gate-method}
\alias{gate<-}
\alias{drawGate}

\title{'gate': a class for subsetting flow-cytometry data by defining
  regions in two-dimesional projections of the data}

\description{In flow-cytometry analysis, regions in two-dimensional
  projections of the data space often have to be selected. Objects of
  this class can store the properties of these selections.}

\section{Creating Objects}{
Objects can be created using methods of the generic function 
\code{\link[prada]{drawGate}} or via\cr
\code{  new("gate",}\cr
\code{    gateFun = ....,   # function returning logical vector}\cr
\code{    colnames = .... # object of class character and length 2}\cr
\code{    logic = ....  # object of class character}\cr
\code{  )}\cr
}

\section{Slots}{\describe{
\item{\code{name}:}{A character vector for the name of the \code{gate}
  object. You can reference the object by its name for subsequent
  operations (e.g. plotting).}
\item{\code{gateFun}:}{A function call together with necessary arguments
  to produce a logical vector when applied on the data.}
\item{\code{colnames}:}{The colnames of the data matrix to which the
gating function is to be applied.}
\item{\code{logic}:}{A \code{character} object, either \code{&} or
  \code{|}. This specifies the logical operation that will be applied
  when combining the selection from the \code{gate} with other object of
that class. See \code{link[prada]{gateSet}} for additional information
on combining \code{gates}.}
\item{\code{type}:}{A \code{character} giving the type of the
  object. This is currently not used but might become important in the
  future.}
}}


\section{Methods}{
  \describe{
    \item{applyGate:}{\code{applyGate(x, data)}\crapplies the gating of
	object \code{x} on data objects of class
      \code{\link[prada]{cytoFrame}} or \code{matrix}.}
    \item{show}{display summary.}
    \item{names}{access slot name.}
}}

\author{Florian Hahne}

\seealso{
  \code{\link[prada]{cytoFrame}}, \code{\link[prada]{gateMatrix}},
\code{\link[prada]{gateSet}}
}

\examples{
g1 <- new("gate", gateFun=function(x)x<500, logic="&",
          colnames=c("FSC-H", "SSC-H"))
}
\keyword{classes}



\name{gateSet-class}
\docType{class}
\alias{gateSet-class}
\alias{gateSet}
\alias{applyGate,gateSet,matrix-method}
\alias{applyGate,gateSet,cytoFrame-method}
\alias{length,gateSet-method}
\alias{show,gateSet-method}
\alias{names,gateSet-method}
\alias{names<-,gateSet-method}
\alias{applyGate}

\title{'gateSet': a class for subsetting flow-cytometry data by defining
  multiple regions in two-dimesional projections of the data}
\description{In flow-cytometry analysis, regions in two-dimensional
  projections of the data space often have to be selected. Objects of
  this class can store the properties for several of these selections}

\section{Creating Objects}{
Objects can be created using methods of the generic function 
\code{\link[prada]{drawGate}} or via\cr
\code{  new("gateSet",}\cr
\code{    glist = ....,   # object of class list}\cr
\code{  )}\cr
}

\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{character} giving the name
    of the object. You can reference the object by its name for subsequent
    operations (e.g. plotting).}
    \item{\code{glist}:}{Object of class \code{"list"} with items of
      class \code{\link[prada]{gate}}. The individual
      \code{\link[prada]{gate}} objects will be combined according to
      the value of their slot \code{logic}.}
  }
}
\section{Methods}{
  \describe{
    \item{applyGate:}{\code{applyGate(x, data)}\crapplies the gating of
      object \code{x} on data objects of class
      \code{\link[prada]{cytoFrame}} or \code{matrix}.}
    \item{length}{length of slot \code{glist}.}
    \item{show}{\code{display summary}
      \item{names}{name of the individual \code{\link[prada]{gate}} objects.}} 
}}

\author{Florian Hahne}

\seealso{
\code{\link[prada]{cytoFrame}}, \code{\link[prada]{gateMatrix}},
\code{\link[prada]{gate}}
}

\examples{
g1 <- new("gate", gateFun=function(x)x<500, logic="&",
          colnames=c("FSC-H", "SSC-H"))
g2 <- new("gate", gateFun=function(x)x>800, logic="&",
          colnames=c("FSC-H", "SSC-H"))
gs <- new("gateSet", name="Set1", glist=list(g1, g2))

}
\keyword{classes}

\name{gateSet-class}
\docType{class}
\alias{gateSet-class}
\alias{gateSet}
\alias{length,gateSet-method}
\alias{show,gateSet-method}
\alias{names,gateSet-method}
\alias{names<-,gateSet-method}
\alias{applyGate}
\alias{[,gateSet,ANY,missing,missing-method}
\alias{[[,gateSet,ANY,missing-method}
\alias{appendGates}
\alias{appendGates, gateSet_method}
\alias{applyGate,matrix,gateSet-method}

\title{'gateSet': a class for subsetting flow-cytometry data by defining
  multiple regions in two-dimesional projections of the data}
\description{In flow-cytometry analysis, regions in two-dimensional
  projections of the data space often have to be selected. Objects of
  this class can store the properties for several of these selections}

\section{Creating Objects}{
Objects can be created using methods of the generic function 
\code{\link{drawGate}} or via\cr
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
      class \code{\link{gate}}. The individual
      \code{\link{gate}} objects will be combined according to
      the value of their slot \code{logic}.}
  }
}
\section{Methods}{
  \describe{
    \item{applyGate:}{\code{applyGate(x, data)} applies the gating of
      object \code{x} on data objects of class
      \code{\link{cytoFrame}} or \code{matrix}}
    \item{length}{length of slot \code{glist}}
    \item{show}{\code{display summary}}
    \item{names, names<-}{extract or replace the names of the
        individual \code{\link{gate}}
	objects.}
    \item{[}{subset to \code{\link{gateSet}}
       objects.}
    \item{[[}{subset to individual \code{\link{gate}}
       objects.}   
    \item{appendGates}{append a \code{gate} or \code{gateSet} to a
       \code{cytoFrame}}
}
}

\author{Florian Hahne}

\seealso{
\code{\link{cytoFrame}}, \code{\link{gateMatrix}},
\code{\link{gate}}
}

\examples{
sampdat <- readFCS(system.file("extdata", "fas-Bcl2-plate323-04-04.A01",
                               package="prada"))
g1 <- new("gate", name="G1", gateFun=function(x)x[,"FSC-H"]<500, logic="&",
          colnames="FSC-H")
g2 <- new("gate", name="G2", gateFun=function(x)x[,"SSC-H"]>800, logic="&",
          colnames="SSC-H")
g3 <- new("gate", name="G3", gateFun=function(x)x[,"FL1-H"]>800, logic="&",
          colnames="FL1-H")
gs <- new("gateSet", name="Set1", glist=list(G1=g1, G2=g2))
length(gs)
gs[[1]]
gs[1]
gsnames <- names(gs)
names(gs) <- gsnames
applyGate(sampdat, gs)
applyGate(exprs(sampdat), gs)
gate(sampdat) <- gs
applyGate(sampdat, 1)
applyGate(sampdat, "G1")
applyGate(sampdat, TRUE)
appendGates(sampdat, g3)

}
\keyword{classes}

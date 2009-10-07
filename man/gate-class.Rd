\name{gate-class}
\docType{class}
\alias{gate-class}
\alias{gate}
\alias{applyGate,matrix,gate-method}
\alias{appendGates,gateSet-method}
\alias{show,gate-method}
\alias{names,gate-method}
\alias{gate<-}
\alias{drawGate}
\alias{combineGates}
\alias{as.gateSet}
\alias{as.gateSet,gate-method}
\alias{names<-,gate-method}
\alias{lines,gate-method}

\title{'gate': a class for subsetting flow-cytometry data by defining
  regions in two-dimesional projections of the data}

\description{In flow-cytometry analysis, regions in two-dimensional
  projections of the data space often have to be selected. Objects of
  this class can store the properties of these selections.}

\section{Creating Objects}{
Objects can be created using methods of the generic function 
\code{\link{drawGate}} or via\cr
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
that class. See \code{link{gateSet}} for additional information
on combining \code{gates}.}
\item{\code{type}:}{A \code{character} giving the type of the
  object. This is currently not used but might become important in the
  future.}  \item{\code{boundaries}:}{A \code{matrix} with two columns
giving the boundaries of the gate in two dimensional space. Can be used
to superimpose the gate boundaries on a plot using \code{lines()}.}
}}


\section{Methods}{
  \describe{
    \item{applyGate:}{\code{applyGate(x, data)} applies the gating of
	object \code{x} on data objects of class
	\code{\link{cytoFrame}} or \code{matrix}.
	In the former case \code{x} may be of class \code{gate},
        \code{gateSet}, \code{character}, \code{numeric} or
        \code{logical}. See vignette for details.}
    \item{show}{display summary.}
    \item{names, names<-}{access and replace slot name.}
    \item{as.gateSet}{Convert \code{gate} object to  \code{gateSet}
      object}
    \item{combineGates}{Combine multiple \code{gate} objects to
       one \code{gateSet} object}
    \item{lines}{Draw the boundaries of the gate.}
 }}

\author{Florian Hahne}

\seealso{
  \code{\link{cytoFrame}}, \code{\link{gateMatrix}},
\code{\link{gateSet}}
}

\examples{
sampdat <- readFCS(system.file("extdata", "fas-Bcl2-plate323-04-04.A01",
                               package="prada"))
g1 <- new("gate", name="test1", gateFun=function(x)x[,"FSC-H"]<500, logic="&",
          colnames="FSC-H", type="misc")
g1
g2 <- new("gate", name="test2", gateFun=function(x)x[,"SSC-H"]>800, logic="&",
          colnames="SSC-H", type="misc")
gs1 <- combineGates(g1,g2)
gs2 <- as.gateSet(g2)
names(g1)
names(g1) <- "testName"
applyGate(sampdat, g1)
applyGate(exprs(sampdat), g2)
gate(sampdat) <- g1
applyGate(sampdat, 1)
applyGate(sampdat, "testName")
applyGate(sampdat, TRUE)
}
\keyword{classes}



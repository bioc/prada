\name{ddCtSet-class}
\docType{class}
\alias{ddCtSet-class}
\alias{exprs,ddCtSet-method}
\alias{exprs<-,ddCtSet,matrix-method}
\alias{initialize,ddCtSet-method}

\title{Class for storage of relative rtPCR results}

\description{A subclass extending virtual class \code{eSet} designed to
  store the raw data and results from a ddCt analysis of rtPCR data.}

\section{Creating Objects}{
  see \code{\link[Biobase]{eSet}} for details
\code{  new("ddCtSet",}\cr
\code{    phenoData = ....,   # object of class \code{AnnotatedDataFrame}}\cr
\code{    experimentData  = .... # object of class \code{MIAME}}\cr
\code{    annotation = ....  # object of class character}\cr
\code{    assayData = ....  # object of class assayData}\cr
\code{  )}\cr
The \code{assayData} slot can be instantiated using function \code{assayDataNew}.
}

\section{Slots}{
  Inherited from \code{eSet}:
   \describe{
      \item{\code{assayData}:}{Contains matrices with equal
	dimensions, and with column number equal to
	\code{nrow(phenoData)}. \code{assayData} must contain matrices
	\code{dCt},\code{dCt.error},\code{Ct.error},
	\code{ddCt},\code{ddCt.error},\code{Difference},\code{level.err},
	\code{number},\code{numberNA},\code{Plate}. See\code{\link[prada]{ddCt}}
	for details.
	Class:\code{\link[Biobase]{AssayData}}}
      \item{\code{phenoData}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{experimentData}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{annotation}:}{See \code{\link[Biobase]{eSet}}}
  }
}

\section{Methods}{
  \describe{
    \item{exprs, exprs<-}{get or set expression data}
    \item{initialize}{constructor, not to be called by the user}
    }}

\author{Florian Hahne}

\keyword{classes}



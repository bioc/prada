\name{cytoFrame-class}
\docType{class}
\alias{cytoFrame-class}

\alias{[,cytoFrame-method}
\alias{[<-,cytoFrame-method}

\alias{description}
\alias{description<-}
\alias{description,cytoFrame-method}
\alias{description<-,cytoFrame-method}

\alias{exprs}
\alias{exprs<-}
\alias{exprs,cytoFrame-method}
\alias{exprs<-,cytoFrame-method}

\alias{show}
\alias{show,cytoFrame-method}
\alias{$.cytoFrame}

\title{Class "cytoFrame": a class for storing observed quantitative
  properties from a population of cells, for example from a FACS run or
  from automated microscopy}

\description{This class represents the data contained in a
  FCS 3.0 file or similar data structures.}

\section{Creating Objects}{
Objects can be created using\cr
\code{  new('cytoFrame,}\cr
\code{    exprs       = ...., # Object of class matrix}\cr
\code{    description = ....  # Object of class character}\cr
\code{  )}\cr
or the function \code{\link{readFCS}}.}


\section{Slots}{
  \describe{
    \item{\code{exprs}:}{Object of class \code{matrix} containing
      the measured intensities. Rows correspond to cells, columns to the
      different channels. The \code{colnames} attribute of the matrix is
      supposed to hold the names or identifiers for the channels. The
      \code{rownames} attribute would usually not be set.}
    \item{\code{description}:}{A named character vector
      containing the experiment description as key-value pairs.}
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "cytoFrame")}: subsetting. Returns
      an object of class \code{cytoFrame}. The subsetting is applied to
      the \code{exprs} slot, while the \code{description} slot is unchanged.}
    \item{[<-}{\code{signature(x = "cytoFrame")}: replaces subsets. }
    \item{exprs}{\code{signature(object = "cytoFrame")}: extract the
      intensities.}
    \item{exprs<-}{\code{signature(object = "cytoFrame")}: replace the
      intensities.}
    \item{description}{\code{signature(object = "cytoFrame")}: extract the
      description.}
    \item{description<-}{\code{signature(object = "cytoFrame")}: replace the
      description.}
    \item{show}{\code{signature(object = "cytoFrame")}: displays a
      summary of the object.}
  }
}

\details{
FCS 3.0 is the Data File Standard for Flow Cytometry, Version FCS3.0.
}
			
\author{Wolfgang Huber}

\seealso{\code{\link{readFCS}}}
\examples{

intens <- matrix(runif(100), ncol=4)
colnames(intens) <- c("FL1-H", "FL2-H", "FL3-H", "FL4-H")

a <- new("cytoFrame",
          exprs=intens,
          description=c(name="example data", date=date()))

description(a)
dim(exprs(a))           

a[1:3, -4]
a[4, 4] <- 4711

}
\keyword{classes}


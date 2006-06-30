\name{readFCS}
\alias{readFCS}
\alias{read.fcs}
\title{Read an FCS 3.0 file}
\description{Read an FCS 3.0 file: Data File Standard for Flow
  Cytometry, Version FCS3.0}
\usage{
  readFCS(filename)
  read.fcs(filename=NULL, objectModel="prada", ...)
}
\arguments{
  \item{filename}{Character of length 1: filename}
  \item{filename}{Character of length 1: the object model to use for the
    output. Either 'prada' for \code{cytoFrame} objects or 'FCS' for
    \code{rflowcyt}'s \code{FCS} objects.
  }
  \item{...}{Arguments that get passed on to funnction \code{read.FCS}}
}

\details{For specifications of FACS 3.0 see
\url{http://www.isac-net.org} and the file
\url{../doc/fcs3.html} in the \code{doc} directory of the package.

The function \code{readFCS} works with the output of the FACS machine
software from a number of vendors. However, the FCS 3.0 standard
includes some options that are not yet implemented in this function. If
you need extensions, please let me know. The output of the fucntion is
an object of class \code{cytoFrame}.

\code{read.fcs} is a wrapper function that allows the user to specify
the class of the output. The purpose of the function is to standardize
the way flow cytometry data is imported into R. Be aware that
\code{rflowcyt} needs to be installed when the function is run with
argument \code{objectModel="FCS"}.  

}

\value{An object of class \code{link[prada]{cytoFrame}} or
  \code{link[rflowcyt]{FCS}}.}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{\code{\link[prada]{readFCSdata}}}

\examples{
sampdat <- readFCS(system.file("extdata", "fas-Bcl2-plate323-04-04.A01",
                               package="prada"))
sampdat
exprs(sampdat[1:3,])
description(sampdat)[3:6]

}
\keyword{IO}


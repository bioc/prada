\name{readFCS}
\alias{readFCS}
\alias{read.fcs}
\title{Read an FCS file}
\description{Read one or several FCS files: Data File Standard for Flow
  Cytometry}
\usage{
  read.fcs(filename=NULL, objectModel="prada", ...)
  readFCS(filename)
}

\arguments{
  \item{filename}{Character of length 1: filename}
  \item{objectModel}{Character of length 1: the object model to use for the
    output. Curently only 'prada' for \code{\link{cytoFrame}} objects is
    supported.}
  \item{...}{Arguments that get passed on to higher-level import functions.}}

\details{
The function \code{readFCS} works with the output of the FACS machine
software from a number of vendors. However, the FCS 3.0 standard
includes some options that are not yet implemented in this function. If
you need extensions, please let me know. The output of the fucntion is
an object of class \code{cytoFrame}.

\code{read.fcs} is a wrapper function that allows the user to specify
the class of the output. The purpose of the function is to standardize
the way flow cytometry data is imported into R using the \code{prada}
package. If the \code{filename} argument to \code{read.fcs} is a
character vector of \code{length > 1}, multiple FCS files can be
imported. Please see the documentation for \code{\link{readCytoSet}} for
alternatives ways to import multiple FCS files and for more details on
the higher-level import function.

For specifications of FCS 3.0 see
\url{http://www.isac-net.org} and the file
\url{../doc/fcs3.html} in the \code{doc} directory of the package.
}

\value{An object of class \code{\link{cytoFrame}}.}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}, Florian Hahne}

\seealso{ \code{\link{readCytoSet}}}

\examples{
sampdat <- readFCS(system.file("extdata", "fas-Bcl2-plate323-04-04.A01",
                               package="prada"))
files <- dir(system.file("extdata", package="prada"),
             pattern="[A-H][0-9][0-9]")
sampdat2 <- read.fcs(system.file("extdata", "fas-Bcl2-plate323-04-04.A01",
                               package="prada"))
sampdat3 <- read.fcs(files, path=system.file("extdata", package="prada"))
sampdat
exprs(sampdat[1:3,])
description(sampdat)[3:6]
class(sampdat3)
}

\keyword{IO}


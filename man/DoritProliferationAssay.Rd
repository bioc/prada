\name{DoritProliferationAssay}

\alias{readWellXML}
\alias{collectdat}
\alias{splitWellName}
\title{Specialized functions for Dorit Arlt's (DKFZ) proliferation 
assay project}
\description{Specialized functions for Dorit Arlt's (DKFZ) proliferation 
assay project}

\usage{
splitWellName(x)
readWellXML(filename, path=".", verbose=TRUE)
collectdat(thedir, verbose=TRUE)
}

\arguments{
  \item{x}{Character of length 1}
  \item{filename}{Character vector of length 1: filename}
  \item{path}{Character vector of length 1}
  \item{verbose}{Logical vector of length 1}
  \item{thedir}{Character.}
}

\details{These functions are rather ideosyncratic for the particular format in
which Urban Liebel's (EMBL) high contents screening microscope stores its data
and how that data was transfered to the DKFZ for Dorit Arlt's proliferation
assay. They are probably not of much interest otherwise.

\code{readWellXML} reads an XML file with the data for one well
from Urban Liebel's (EMBL) HCS-microscope.

\code{splitWellName} splits a file name string into its components.

\code{collectdat} collects together data from the individual wells 
into one big data frame.
}

\value{Various.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\seealso{}
\examples{
  splitWellName("dorit10_3c_1-2--01-pDhfbr2_78k24_1E3Y1")
  sampdat <- readWellXML(file="dorit13_2c-1-2--13-Ytrail=3143937171.xml",
                   path=system.file("extdata", package="prada"))
  sampdat
}
\keyword{IO}


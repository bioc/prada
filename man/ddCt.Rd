\name{ddCt}
\alias{ddCt}
\title{Apply the ddCt method for a given data set}
\description{Apply the ddCt method for a given data set}
\usage{
ddCt(raw.table,calibrationSample,housekeepingGene,sampleInformation=NULL,filename="warning.output.txt")
}
\arguments{
  \item{raw.table}{data frame. It must contain columns with the
    following names:'Ct','Sample','Detector','Platename'. The column
    'Ct' must contain numeric values.}
  \item{calibrationSample}{character. Name of the calibration samples.}
  \item{housekeepingGene}{character of length 1: name of the
    housekeeping gene.}
  \item{sampleInformation}{if specified it must be an object of class
    \code{phenoData} with a column named 'Sample'.}
  \item{filename}{character of length 1. The name of the file the
    warnings should be stored in.}
}

\details{Please see vignette 'How to apply the ddCt method' in the prada
  package for details.}

\value{A an object of class
  \code{\link[Bioconductor:eSet-class]{eSet}}.The names of the eList
  have the following names: exprs(includes the level),
  Ct,dCt,ddCt,MAD,Difference,number\_NA,Plate.}

\seealso{\code{\link{readSDM}},\code{\link{eListWrite}}}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\keyword{}
\examples{##see vignette rtPCR.pdf}

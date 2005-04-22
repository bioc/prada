\name{ddCt}
\alias{ddCt}
\title{Apply the ddCt method for a given data set}
\description{Apply the ddCt method for a given data set}
\usage{
ddCt(raw.table,calibrationSample,housekeepingGenes,type="mean",sampleInformation=NULL,toZero=FALSE,filename="warning.output.txt")
}

\arguments{
  \item{raw.table}{data frame. It must contain columns with the
    following names:'Ct','Sample','Detector','Platename'. The column
    'Ct' must contain numeric values.}
  \item{calibrationSample}{character. Name of the calibration samples.}
  \item{housekeepingGenes}{character. Name of the housekeeping genes.}
  \item{sampleInformation}{if specified it must be an object of class
    \code{phenoData} with a column named 'Sample'.}
  \item{filename}{character of length 1. The name of the file the
    warnings should be stored in.}
  \item{type}{character of length 1. `mean` or `median`- which method should be used for the
    aggregation of the repicates}
   \item{toZero}{boolean - if there is only one replication should the error be
    treated as zero ? (only if 'type' is mean)}
}

\details{Please see vignette 'How to apply the ddCt method' in the ddCt
  package for details.}

\value{A an object of class
  \code{\link[Bioconductor:eSet-class]{eSet}}. Its eList
  has the following names: exprs(the level), level.err,
  Ct, Ct.error, dCt, dCt.error, ddCt, ddCt.error, Difference,
  number\_NA, number, Plate.}

\seealso{\code{\link{readSDM}},\code{\link{eListWrite}}}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\keyword{}
\examples{##see vignette rtPCR.pdf}

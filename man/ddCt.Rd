\name{ddCt}
\alias{ddCt}
\title{Apply the ddCt method for a given data set (with or without efficiencies}
\description{Apply the ddCt method for a given data set}
\usage{
ddCt(raw.table,calibrationSample,housekeepingGenes,type="mean",sampleInformation=NULL,filename="warning.output.txt")
}

\arguments{
  \item{raw.table}{data frame. It must contain columns with the
    following names:'Ct','Sample','Detector','Platename'. The column
    'Ct' must contain numeric values.}
  \item{calibrationSample}{character. Name of the calibration samples.}
  \item{housekeepingGenes}{character of length 1: name of the
    housekeeping gene.}
  \item{efficiencies}{named vector with efficiencies for all genes}
  \item{efficiencies.error}{named vector with error for efficiencies for
    all genes}
  \item{sampleInformation}{if specified it must be an object of class
    \code{phenoData} with a column named 'Sample'.}
  \item{filename}{character of length 1. The name of the file the
    warnings should be stored in.}
  \item{type}{`mean` or `median`- which method should be used for the
    aggregation of the repicates}
}

\details{Please see vignette 'How to apply the ddCt method' in the ddCt
  package for details.}

\value{A an object of class
  \code{\link[Bioconductor:eSet-class]{eSet}}.For ddCt the eList
  has the following names: exprs(the level),level.err,
  Ct,Ct.error,dCt,dCt.error,ddCt,ddCt.error,Difference,number\_NA,number,Plate.
  For ddCtwithE the eList has the following names: exprs(the level),level.err,
  Ct,Ct.error,efficiencies,efficiencies.error,Difference,number\_NA,number,Plate.}

\seealso{\code{\link{readSDM}},\code{\link{eListWrite}}}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\keyword{}
\examples{##see vignette rtPCR.pdf}

\name{statWell}
\alias{statWellLocfit}
\title{Calculate a statistic for the data from one well}
\description{Calculate a statistic for the data from one well}
\usage{
statWellLocfit(x, span,
  plotwhat="nothing", plotdir=".", plotfile, ...)
}
\arguments{
  \item{x}{data.frame with (at least) columns \code{brdu}, \code{trsf},
    \code{dapi}, \code{Field}, and \code{cloneId}.}
  \item{span}{numeric. Is used as the parameter of \code{alpha}
    in the call to the function \code{\link[locfit:locfit]{locfit}}.}
  \item{plotwhat}{character of length 1. Allowed values are
    \code{nothing}, \code{screen}, \code{figscp}, \code{mkpp}.
    The parameter is used to control the number and format of the
    diagnostic plots that are produced as a side effect of
    this function.}
  \item{plotdir}{character of length 1. The name of directory where
    diagnostic plots will be saved (only relevant if \code{plotwhat} is
    \code{figscp} or \code{mkpp}).}
  \item{plotfile}{character. The file name root for the file with
    the diagnostic plot (only relevant if \code{plotwhat} is
    \code{figscp} or \code{mkpp}). If missing, a default name is
    generated.}
  \item{...}{Further arguments that are passed on to \code{\link{plot}}.}
}

\details{The function fits a local regression model to the data. For
  details, please see the paper.}

\value{A list with elements \code{cloneId}, \code{nrcells}, \code{trsfeff},
  \code{delta}, \code{se.delta}, \code{zscore}, \code{plotfile}.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}
\examples{##see vignette}

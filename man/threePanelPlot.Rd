\name{threePanelPlot}
\alias{threePanelPlot}
\title{Visualize cytometry data}
\description{
  Function to visualize multivariate (cytometry) data in three
  two-dimensional plots.
}
\usage{
threePanelPlot(data, x.panels = c(1, 4, 5), y.panels = c(2, 3, 6),
               tot.width = 15, tot.height = 5.4, maxcells = 20000,
               limits = c(0, 1023), remove.extremes = TRUE,
               plotTitle = "Three-Panel Plot", use.smoothScatter = TRUE,
               palette = colorRampPalette(brewer.pal(9, "Blues")),
               new.device = TRUE, verbose = TRUE,
               addPoints = NULL, addCol = "red", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data matrix to visualize}
  \item{x.panels}{which variables (columns) are to be plotted at the
    x-axis of the three variables}
  \item{y.panels}{which variables (columns) are to be plotted at the
    y-axis of the three variables}
  \item{tot.width}{width of a new device to open, see argument
    \code{new.device}}
  \item{tot.height}{height of a new device to open, see argument
    \code{new.device}}
  \item{maxcells}{maximum number of observations (cells) for plotting;
    higher numbers reduce performance}
  \item{limits}{minimum and maximum value (theoretically) observed in
    the data; e.g., with 10-channel digitized data it is c(0,1023)}
  \item{remove.extremes}{logical; are extreme values (equal to
    theoretical \code{limits}) to be removed before plotting}
  \item{plotTitle}{title for the plot}
  \item{use.smoothScatter}{logical, should the function
    \code{\link[graphics]{smoothScatter}} be employed for plotting the data
    (plots data densities rather than individual points)}
  \item{palette}{if \code{smoothScatter} is used, which colour palette
    is it to use}
  \item{new.device}{logical; should a new device be opened for the three
    plots; if \code{FALSE} the three plots will be plotted to the
    currently active device}
  \item{verbose}{logical; do you want extended output to STDOUT}
  \item{addPoints}{should special points be marked after plotting the
    data; is expected to be a subser of argument \code{data} with the
    same number of columns (=variables); if \code{NULL} no points are marked}
  \item{addCol}{in which colour are the points in \code{addPoints} to be
  marked}
  \item{\dots}{further arguments passed on to \code{plot.default}}
}
%\details{~~ If necessary, more details than the __description__  above ~~}
\value{
  no value is returned; the function is called to produce three plots
}
%\references{ ~put references to the literature/web site here ~ }
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
%\note{ ~~further notes~~ }
\seealso{\code{\link{plot.default}}}
\examples{
   # generate some data:
   toyData <- cbind(matrix(pmax(0,pmin(runif(3000)+rnorm(3000),4)),ncol=3),
                    matrix(pmax(0,pmin(rnorm(3000,2,1),4)),ncol=3))
   colnames(toyData) <- paste("Var",1:6,sep="")
   toyQuantiles <- apply(toyData,2,quantile,probs=c(0.25,0.5,0.75))

   # plot it and mark the quantiles:
   threePanelPlot(toyData,addPoints=toyQuantiles,
                  addCol=c("orange","red","purple"),limits=c(0,4),pch=20)
}
\keyword{hplot}% at least one, from doc/KEYWORDS


\name{plotGate}
\alias{plot.cytoFrameGate}
\alias{plot.cytoSetGate}
\alias{plot.matrixGate}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualize Observations within defined Gates}
\description{
  This functions visualize which observations are included in gates
  defined by the functions \code{gate.cytoFrame}, \code{gate.cytoSet},
  and \code{gate.matrix}. The make use of the visualization function
  \code{\link[prada]{threePanelPlot}} .
}
\usage{
plot.cytoFrameGate(x, y, x.panels = c(1, 4, 5), y.panels = c(2, 3, 6),
                   use.smoothScatter = TRUE, limits = c(0, 1023),
                   plotTitle = "Gated Data", gate.colour = "red",
                   new.device = TRUE, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{gate resulting, e.g.,  from function
    \code{\link[prada]{gate.cytoFrame}}; thus \code{x} is of class
      \code{cytoFrameGate} or similar}
  \item{y}{the original data used to derive gate \code{x}, thus either a
    \code{cytoFrame}, \code{cytoSet} or \code{matrix}}
  \item{x.panels}{which variables (columns) are to be plotted at the
    x-axis of the three variables}
  \item{y.panels}{which variables (columns) are to be plotted at the
    y-axis of the three variables}
  \item{use.smoothScatter}{logical, should the function
    \code{\link[geneplotter]{smoothScatter}} be employed for plotting the data
    (plots data densities rather than individual points)}
  \item{limits}{minimum and maximum value (theoretically) observed in
    the data; e.g., with 10-channel digitized data it is c(0,1023)}
  \item{plotTitle}{title for the plot}
  \item{gate.colour}{which colour to use for observations within the
    gate \code{x}}
  \item{new.device}{logical; should a new device be opened for the three
    plots; if \code{FALSE} the three plots will be plotted to the
    currently active device}
  \item{verbose}{logical; do you want extended output to STDOUT}
  \item{\dots}{further arguments passed on to \code{plot.default}}
}
%\details{}
\value{
  No value is returned. These functions are called for producing the
  three plots.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
%\note{ ~~further notes~~ }
\seealso{\code{\link[prada]{gate.cytoFrame}},\code{\link[prada]{threePanelPlot}}}
\examples{
   # see the example in function 'gate.cytoFrame'
}
\keyword{hplot}% at least one, from doc/KEYWORDS


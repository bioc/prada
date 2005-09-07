\name{gate}
\alias{gate.cytoFrame}
\alias{gate.matrix}
\alias{gate.cytoSet}

\title{Functions to gate cytometry data}
\description{
  These functions allow the user to select a number of observations in
  multivariate data by circling them in two-dimensional projections of
  the data. \code{gate.cytoSet} calls \code{gate.cytoFrame} seperately
  on each frame in the set. \code{gate.matrix} does the same as
  \code{gate.cytoSet} but only requires a data matrix.
}
\usage{
my.gate <- gate.cytoFrame(object, gate.colour = "red",
                          use.smoothScatter = FALSE,
                          data.min = 0, data.max = 1023,
                          max.observations = 20000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{cytoFrame}, \code{cytoSet} or matrix of cytometry
    data, from which observations are to be selected by gating}
  \item{gate.colour}{colour of the drawn polygons}
  \item{use.smoothScatter}{logical, should the function
    \code{\link{smoothScatter}} be employed for plotting the data
    (plots data densities rather than individual points)}
  \item{data.min}{minimum value to (theoretically occur in the data;
    e.g. with 10-channel-digitized data, it is 0}
  \item{data.max}{minimum value to (theoretically occur in the data;
    e.g. with 10-channel-digitized data, it is 1023}
  \item{max.observations}{maximum number of observations from the data
    to be plotted; large values decrease performance}
}
\details{
  The user is prompted to select the two dimensions
  for the projection and the defines a polygon by left-mouse clicks on
  its vertices. To finish drawing the polygon, press the right mouse
  button. Observations within the polygon are computed and marked. The
  user is then prompted if he wants to
  \itemize{
    \item redraw the polygon
    \item combine this polygon with gated region to generate
      unions of regions
    \item proceed, in which case only observations within the selected
      region(s) are kept, while the rest is discarded. The thus selected
      observations can then be further reduced by another gating in
      other two dimensions OR
    \item finish, in which case observations within the final gate
      are marked as TRUE, while the rest is set to FALSE
  }
}
\value{
  For \code{gate.cytoFrame} and \code{gate.matrix}, it is a logical
  vector stating for each observation in the data (in same order) wether
  it was included in the defined gate (\code{TRUE}) or not
  (\code{FALSE}). For \code{gate.cytoSet}, it is a list of logical
  vectors, one for each \code{cytoFrame} in the \code{cytoSet}.
  The class of the result depends on the called function and is
  \code{cytoFrameGate},\code{cytoSetGate}, or \code{matrixGate}. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{Classes \code{\link{cytoFrame-class}},
  \code{\link{cytoSet-class}}, \code{\link{plot.cytoFrameGate}} for
  plotting}
\examples{
  \dontrun{
    # don't run during package test, because requires user input

    # read in FCS files:
    mySet   <- readCytoSet(...)

    # let user do gating:
    mySGate <- gate.cytoSet(mySet)

    # visualize gating:
    plot(mySGate,mySet,pch=".",cex=2)
  }
}
\keyword{iplot}% at least one, from doc/KEYWORDS

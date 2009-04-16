\name{plotPlate}
\alias{plotPlate}
\title{Plot a well statistic for microtiter plates.}
\description{Plot a well statistic in false color representation or
  using a self-defined grid plotting function. The
  plot is supposed to resemble the physical geometry of a microtitre
  plate.}

\usage{
plotPlate(x,nrow = 8, ncol = 12, col=c("red", "blue"),
ind = 1:(ncol*nrow), xrange=range(x, na.rm=TRUE), na.action = "zero",
main, char, desc = character(2), add=FALSE, gridFun="default",
funArgs=NULL,...)
}

\arguments{
  \item{x}{Numeric vector of length \code{ncol*nrow} or matrix with
    \code{ncol*nrow} rows (except if argument \code{ind} is specified).
    If of class \code{matrix}, the use of argument \code{gridFun} is expected.}

  \item{nrow}{Numeric of length 1. The number of rows of the plate.}

  \item{ncol}{Numeric of length 1. The number of columns of the plate.}

  \item{col}{Character vector. Usually the names of two or three colors
    between which the color map is interpolated, using the function
    \code{\link[grDevices::colorRampPalette]{colorRampPalette}}.}

  \item{ind}{Optional integer vector of equal length as \code{x}. It
    indicates the position of the respective value of
    x on the plate. Can be used to adress the problem of missing
    values. Each well that is not allocated a value of \code{x} by
    \code{ind} will not be plotted.}

  \item{xrange}{Range of \code{x} that is mapped into the color scale.}

  \item{na.action}{Character. One of \code{"zero"} \code{"omit"} or
    \code{"xout"}.
    How should the wells for which \code{x} is NA be treated?
    For \code{"zero"}, they are plotted as if the value were 0. 
    For \code{"omit"}, they are omitted.
    For \code{"xout"}, they are crossed out.
    When \code{x} is a matrix, \code{na.action} is only applied to rows
    containing nothing but NAs. Further special treatment of NA
    values in matrices need to be implemented in \code{gridFun}.}

  \item{main}{Character of length 1. Plot title.}

  \item{char}{An optional character vector of equal length as \code{x}
    (except if argument \code{ind} is specified) to be used for well
    annotation. Each element of the
    vector may contain a string to be superimposed on
    the respective well or \code{NA} for no plotting.}
 
  \item{desc}{Character of length 2. Legend for the two extremes of the
    colorbar, e.g. 'act' and 'inh'. }
 
  \item{add}{Logical. If \code{TRUE} add plate plot to current plot. May
    be used when plotting in grid layout panels.} 
    
  \item{gridFun}{Character. The name of the plotting function to create
    individual graphs for each well. See functions \code{.drawCircle}
    and \code{.drawPie}} for examples.

  \item{funArgs}{Dataframe with argument values to be passed to
    gridCall. For each argument specified in gridCall there must be one column
    with the argument name as colname and the argument values for every
    well.}

  \item{...}{Further graphical parameters that can be used to control
    the output of plotPlate.
    \describe{
    \item{cex.main:}{expansion factor for title.}
    \item{cex.lab:}{expansion factor for label.}
    \item{cex.char:}{expansion factor for well annotation.}
    \item{cex.legend:}{expansion factor for well legend labels.}
    \item{cex.desc:}{expansion factor for well legend description.}
  }
  }
}

\details{ You may use this function either to create plots showing a
single-value per well statistic for microtiter plates, or you can use a
self-made plotting function using a combination of any valid grid
commands to produce arbitrary plots in a plate array format. These plots
may also show multifactorial data. Self-defined plotting functions need
to have \code{data} as first argument. \code{plotPlate} passes all data
values for the respective well to the plotting function. Any further
arguments may be passed on using argument \code{funArgs}. See
\code{.drawCircle} and \code{.drawPie} for examples of valid plotting
functions and the vignette for detailed information.  Note that
using \code{funCall} overrides some of the default functionalities,
e.g. plotting of legends and alters the treatment of NA values.

Argument \code{ind} allows the user to indicate the position (well
number) for each element of vector \code{x} on the plate. This can be
used either to change the order in which elements of \code{x} are to be
plotted or to deal with the problem of missing data for some of the
wells on a plate.

To further increase the amount of information of the platePlot one may
decorate wells with short annotations using argument \code{char}. Each
element of \code{char} \code{!= NA} will be superimposed on the respective
well (see examples). 
}

\value{The function produces a plot in the active graphics device.  

It returns a list withfour elements. The element \code{which} is a vector with
the indices of those elements in \code{x} that were plotted (see argument
\code{na.action}). The element \code{coord} is a \code{length(which)} by 4
matrix in which each row specifies the corners of a rectangle that contains a
well. It is intended to be used as an argument to a subsequent call to
\code{\link[geneplotter::imageMap]{imageMap}}. Elements \code{width} and
\code{height} may be used to open a graphic devices that can hold the plate
plot with the correct aspect ratio.}

\seealso{\code{\link[geneplotter::imageMap]{imageMap}}}

\author{Florian Hahne, Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\keyword{hplot}

\examples{
plotPlate(runif(96), main="example 1", col=c("#0000e0", "#e00000"), desc=c("act", "inh"))
plotPlate(runif(384), nrow=16, ncol=24, main="example 2", col=c("#0000e0", "white", "#e00000"))
plotPlate(runif(48), main="example 3", col=c("#0000e0", "#e00000"), ind=c(1:24, 73:96))
x <- runif(96)
x[sample(96, 10)] <- NA
plotPlate(x, main="example 4", col=c("#0000e0", "#e00000"),
char=c(rep(NA, 72), LETTERS[1:24]), na.action="xout")
plotPlate(runif(96, min=0.1, max=0.5), gridFun=".drawCircle")
plotPlate(matrix(runif(288), ncol=3), gridFun=".drawPie",
funArgs=as.data.frame(matrix(2:4, ncol=3, nrow=96, byrow=TRUE)))
 
}

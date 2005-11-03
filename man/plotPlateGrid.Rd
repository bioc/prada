\name{plotPlateGrid}
\alias{plotPlateGrid}
\title{Plot a well statistic in false color representation.}
\description{Plot a well statistic using any tailored grid plotting function. The
  plot is supposed to resemble the physical geometry of a microtitre
  plate.}

\usage{
plotPlateGrid(x, gridCall="circle", callArgs=NULL, nrow = 8, ncol = 12, ind = 1:(ncol*nrow), main, xrange,
          col, device="x11", file, width, na.action = "zero",
          desc = as.character(c(NA, NA)), char, legend = TRUE)
}

\arguments{
  \item{x}{Numeric vector of length \code{ncol*nrow} or matrix with
    \code{nrow=ncol*nrow} (except if argument
    \code{ind} is specified).}
  \item{gridCall}{Character. The name of the plotting function to create
    individual graphs for each well.}
  \item{callArgs}{Dataframe with argument values to be passed to
    gridCall. For each argument specified in gridCall you need a column
    with the argument name as colname and the argument values for each well.}
  \item{nrow}{Numeric of length 1. The number of rows of the plate.}
  \item{ncol}{Numeric of length 1. The number of columns of the plate.}
  \item{ind}{Optional integer vector of equal length as \code{x}. It
    indicates the position of the respective value of
    x on the plate. Can be used to adress the problem of missing
    values. Each well that is not allocated a value of \code{x} by
    \code{ind} will not be plotted.}
  \item{main}{Character of length 1. Plot title. This arguments gets passed 
     on to \code{\link[graphics:text]{text}} as the argument \code{labels}.}
  \item{xrange}{Range of \code{x} that is mapped into the color scale.}
  \item{col}{Character vector. Usually the names of two or three colors
    between which the color map is interpolated, using the function
    \code{\link[geneplotter::colorRampPalette]{colorRampPalette}}.}
  \item{device}{Character. Name of a valid graphics device,
    e.g. \code{x11}, \code{png}, \code{pdf}.}
  \item{file}{The \code{file} argument of the call to
    \code{device}, i.e. the name of the graphics output file.}
  \item{width}{The \code{width} argument of the call to 
    \code{device}, i.e. the width of the graphic.}
  \item{na.action}{Character. One of \code{"zero"} \code{"omit"} or
    \code{"xout"}.
    How should the wells for which \code{x} is NA be treated?
    For \code{"zero"}, they are plotted as if the value were 0. 
    For \code{"omit"}, they are omitted.
    For \code{"xout"}, they are crossed out.}
  \item{desc}{Character of length 2. Legend for the two ends of the
    colorbar, e.g. 'act' and 'inh'. }
  \item{char}{An optional character vector of equal length as \code{x}
    (except if argument \code{ind} is specified). Each element of the
    vector may contain one single \code{char} to be superimposed on
    the respective well or \code{NA} for no plotting.}
  \item{legend}{Logical. Should legend be plotted?}
}

\details{
You may use this function either to create the same plots like in 
\code{\link[prada::plotPlate]{plotPlate}}, or you can use a self-made
plotting function using a combination of any valid grid commands. This
function needs to have \code{data} as first argument. Any further
arguments may be passed on using \code{callArgs}. See 
\code{.drawCircle} and \code{.drawPie} for examples of valid plotting
functions. Using gridCall overrides some of the default fucntionalities,
e.g. plotting of a legend.

Argument\code{ind} allows the user to indicate the position (well
number) for each element of vector \code{x} on the plate. This can be
used either to change the order in which elements of \code{x} are to be
plotted or to deal with the problem of missing data for some of the
wells on a plate.

To further increase the amount of information of the platePlot one may
decorate wells with single characters using argument \code{char}. Each
element of \code{char} \code{!= NA} will be superimposed on the respective
well (see examples). 
}

\value{The function produces a plot in the graphics device specified by the
\code{device} argument.  

It returns a list with two elements. The element \code{which} is a vector with
the indices of those elements in \code{x} that were plotted (see argument
\code{na.action}). The element \code{coord} is a \code{length(which)} by 4
matrix in which each row specifies the corners of a rectangle that contains a
well. It is intended to be use as an argument to a subsequent call to
\code{\link[geneplotter::imageMap]{imageMap}}.}

\seealso{\code{\link[prada::imageMap]{imageMap}}}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\keyword{hplot}

\examples{
if((capabilities("X11") || options("device")=="windows") &
interactive()){

x <- matrix(runif(96*3), ncol=3)
y <- runif(96)

plotPlateGrid(y, width=7, col=5:6, main="exampe 1")

plotPlateGrid(y, ".drawCircle", width=7, main="exampe 2")

callArgs <- data.frame(col1=rep(2,96), col2=rep(3,96), col3=rep(4,96))
plotPlateGrid(x, ".drawPie", callArgs, width=7, main="exampe 3")

}
}

\name{imageMap}
\alias{imageMap}
\title{Write an HTML IMG tag together with a MAP image map.}
\description{Write an HTML IMG tag together with a MAP image map.}
\usage{
  imageMap(con, imgname, coord, tooltips, url, target="extra")
}
\arguments{
  \item{con}{Connection (see argument \code{con} of \code{\link{writeLines}}).}
  \item{imgname}{Character. Name of the image file for which you want to
    construct the map.}
  \item{coord}{Matrix with 4 columns. Each row specifies the corners of a 
  rectangle within the image, in the following order: (left x,
  lower y, right x, upper y).}
\item{tooltips}{Character vector of length \code{nrow(coord)}. It
  contains the text that is displayed when the mouse moves over a rectangle.}
\item{url}{Character vector of length \code{nrow(coord)}. It contains
  the hyperlinks that the browser can follow when the mouse clicks into
  a rectangle.}
  \item{target}{Character. Name of the target browser window in which
    the target of the hyperlinks are displayed. If its length is less
    than \code{nrow(coord)}, the recycling rule is invoked.}
}

\details{This function is typically used within the following sequence
of steps:
\enumerate{
\item generate your graphic and save it as a bitmap file, e.g.
  using the \code{\link{jpeg}} or \code{\link{png}} devices. At this
  stage, you also need to figure out the pixel coordinates of the
  interesting regions within your graphic. Since the mapping between
  device coordinates and pixel coordinates is not very transparent,
  this may be a little tricky. See the examples below, and for a more
  complex example, see the source code of the function
  \code{\link{plotPlate}}.
\item open an HTML page for writing, e.g. using the
  \code{\link{openHTMLpage}} function.
\item optionally write some header information into the HTML page.
\item call the \code{\link{imageMap}} function, using for \code{coord}
  the coordinates you have from step 1.
\item optionally write some further information into the HTML page.
\item close the HTML page, e.g. using the
  \code{\link{closeHTMLpage}} function.
  }
}

\value{The function is called for its side effect, which is writing text into
the connection \code{con}.}

\seealso{\code{\link[prada::plotPlate]{plotPlate}}, 
   \code{\link[base::writeLines]{writeLines}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}
\examples{
fhtml  <- tempfile()
fpng   <- tempfile()

if(capabilities()["png"]) {
  width <- height <- 512
  png(fpng, width=width, height=height)
  par(mai=rep(0,4))
  plot(0, xlim=c(0, width-1), ylim=c(0, height-1), xaxs="i", yaxs="i",
       type="n")
  coord=cbind(c(100,360),c(400, 50),c(200,460),c(500,150))
  rect(coord[,1], height-coord[,2], coord[,3], height-coord[,4], col=c("#1f77b4", "#b2df8a"))
  dev.off()

  con <- openHTMLpage(fhtml, "Hello world")
  imageMap(con, fpng, coord=coord, 
    url=paste("http://", c("link1", "link2"), sep=""),
    tooltips=c("Blue box", "Green box"))
  closeHTMLpage(con)

  cat("Now have a look at file ", fhtml, ".html with your browser.\n", sep="")
}
}



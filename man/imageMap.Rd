\name{imageMap}
\alias{imageMap}
\title{Write an HTML IMG tag together with a MAP image map.}
\description{Write an HTML IMG tag together with a MAP image map.}
\usage{
  imageMap(con, imgname, coord, tooltips, url, target="extra")
}
\arguments{
  \item{con}{Connection (see argument \code{con} of \code{\link{writeLines}}).}
  \item{imgname}{Character. Name of the image file, as it is to appear in the
  HTML output.}
  \item{coord}{Matrix with 4 columns. Each row specifies the corners of a 
  rectangle within the image.}
  \item{tooltips}{Character of length \code{nrow(coord)}.}
  \item{url}{Character of length \code{nrow(coord)}.}
  \item{target}{Character. Name of the target browser window.}
}

\details{See example.}

\value{The function is called for its side effect, which is writing text into
the connection \code{con}.}

\seealso{\code{\link[prada::plotPlate]{plotPlate}}, 
   \code{\link[base::writeLines]{writeLines}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}
\examples{
imageMap(stdout(), "myimage.jpg", coord=matrix(1:8,nrow=2), 
   url=c("a","b"), tooltips=c("TT1", "TT2"))
}



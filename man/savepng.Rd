\name{savepng}
\alias{savepdf}
\alias{savepng}
\title{Save the contents of the current graphics device to a PDF or PNG file}
\description{Save the contents of the current graphics device to a PDF or PNG file}

\usage{
  savepng(fn, dir, width=480, asp=1)
  savepdf(fn, dir, width=6, asp=1)
}

\arguments{
  \item{fn}{character: name of the output file (without extension .png or .pdf)}
  \item{dir}{character: directory to which the file should be written.}
  \item{width}{numeric: width of the image in pixels (png) or inches (pdf)}
  \item{asp}{numeric: aspect ratio; height=width*asp}
}

\details{}

\value{Character: name of the file that was written.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}

\seealso{\code{\link{dev.copy}},\code{\link{pdf}},\code{\link{png}}}

\examples{
  x = seq(0, 20*pi, len=1000)
  plot(x*sin(x), x*cos(x), type="l")

  try({   ## on some machines, png or pdf may not be available
  savepdf("spiral", dir=tempdir())
  savepng("spiral", dir=tempdir())
  })    
}

\keyword{programming}
\keyword{error}

\name{openHTMLpage}
\alias{openHTMLpage}
\alias{closeHTMLpage}
\title{Open and close an HTML file for writing.}
\description{Open and close an HTML file for writing..}
\usage{
  openHTMLpage(name, title="")
  closeHTMLpage(con)
}
\arguments{
  \item{name}{Character. File name (\emph{without} the extension '.html').}
  \item{title}{Character. Value of the \code{title} tag in the HTML header.}
  \item{con}{Connection.}
}

\details{See example.}

\value{For \code{openHTMLpage}, a \code{\link[base::connection]{connection}}. }

\seealso{}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{IO}
\examples{
  fn <- tempfile()
  con <- openHTMLpage(fn, "My page")
  writeLines("Hello world", con)
  closeHTMLpage(con)
  readLines(paste(fn, ".html", sep=""))
}



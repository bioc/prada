\name{getPradaPar}
\alias{setPradaPars}
\alias{getPradaPar}
\title{Set and query global parameters for functions in the prada package}
\description{Set and query global parameters for functions in the prada package}
\usage{
setPradaPars(pars)
getPradaPar(parname)
}
\arguments{
\item{pars}{Named list}
\item{parname}{Character string of length 1}
}
\details{TBA}

\value{For \code{getPradaPar}, the value of the list element with name
\code{parname} of the global parameters list.
The function \code{setPradaPars} is invoked for its side effect,
which is assigning a value to the global parameters list.
It returns the value invisible(NULL).
}
\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\examples{
## setPradaPars(list(tomato=1, apple="two", pear=as.integer(3)))
## getPradaPar("pear")
}
\keyword{manip}


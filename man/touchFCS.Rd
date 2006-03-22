\name{touchFCS}

\alias{touchFCS}

\title{Check for FCS files}

\description{
  The function reads the header of a file or of a range of files and
  checks wether they are valid FCS 2.0 or FCS 3.0 files.
}

\usage{
touchFCS(path = ".", file)
}

\arguments{
  \item{path}{character, the path to a folder containing files}
  \item{file}{character, the path to a single file}
}

\details{
  The user may either specify the path to a directory in which to search
  for FCS files or the path to a single file. 
}

\value{
 A character vector with names of the valid FCS files found. 
}

\author{fhahne}



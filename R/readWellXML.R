readWellXML = function(file, path=".", verbose=TRUE) {
  stopifnot(is.character(file) && length(file)==1)
  stopifnot(is.character(path) && length(path)==1)
  stopifnot(is.logical(verbose) && length(verbose)==1)
  
  if(verbose) cat(file, "\t")
  z = xmlTreeParse(file.path(path, file))
  r = xmlRoot(z)[[2]]
  ## Name, NumElts, String, Array 
  stopifnot(xmlValue(r[[2]]) == "2")
  
  ## wellname
  wellname = getElement(r[[3]], name="current clone", type="character")

  ## check whether consistent with file
  f = strsplit(file, "=")[[1]][1]
  stopifnot(f==wellname)  
  welldat = splitWellName(f)
  if(verbose) cat(welldat, sep="\t")
 
  ## the actual data array
  ch = xmlChildren(r[[4]])
  stopifnot(xmlValue(ch[[1]]) == "reports cluster")

  ## the data come as a matrix (with dimensions dims) of vectors
  ## dims[1] is the number of microscope fields (images)
  ## dims[2] is the maximum number of cells in a microscope field
  ## In those fields that do not contain the maximum number of cells
  ## (which are most), the values are padded with zeros.
  for (i in 2:3) stopifnot(xmlName(ch[[i]]) == "Dimsize")
  dims = as.numeric(c(xmlValue(ch[[2]]), xmlValue(ch[[3]])))
  stopifnot(length(ch)-3 == prod(dims))

  thenames = list("Area (pixels)", "Perimeter", "Projection x", "Projection y",
      "Longest Segment Length",
      list("x1Left", "y1Top", "x2Right", "y2Bottom"),
      "Number of Holes", "Hole's Area (pixels)", "Dapi mean",
      "Brdu mean", "transfection mean")
  thetypes = c(rep("numeric", 5), list(rep("numeric", 4)), rep("numeric", 5))

  znames      = unlist(thenames)
  y           = matrix(NA, nrow=prod(dims), ncol=length(znames)+1)
  colnames(y) = c(znames, "Field")
  for(i in 1:nrow(y)) {
    if(verbose && (i%%50 == 0)) cat(".")
    z = getVector(ch[[3+i]], name="Cluster", elementnames=thenames, elementtypes=thetypes)
    stopifnot(all(names(z)==znames))
    y[i, znames] = z
  }

  ## remove dummy columns
  dummylines = (rowSums(abs(y), na.rm=TRUE) == 0) 
  y[, "Field"] = 1 + ((0:(nrow(y)-1)) %/% dims[2])
  ny = matrix(y[!dummylines, ], ncol=ncol(y))
  colnames(ny) = colnames(y)

  if(verbose) cat("", nrow(ny), "cells\n")

  wdrep = matrix(welldat, ncol=length(welldat), nrow=nrow(ny), byrow=TRUE)
  colnames(wdrep) = names(welldat)
  
  return(cbind.data.frame(wdrep, ny))
}

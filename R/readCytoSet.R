readCytoSet <- function(files, path=".", ...) {
  ## if(!is.(pd, "phenoData"))
  ##  stop("'pd' must be an object of class 'phenoData'")
  if(missing(files)) {
    files <- dir(path, ...)
    if(length(files)<1)
      stop(paste("No (matching) files found in directory", path))
  }
  if(!is.character(files))
    stop("'files' must be a character vector")

  env <- new.env(hash=TRUE)
  cn  <- NULL
  for(f in files) {
    x <- readFCS(file.path(path, f))
    if(is.null(cn)) {
      cn <-colnames(x)
    } else {
      if(!identical(cn, colnames(x)))
        stop(paste(f, "has different 'colnames' than previously read files."))
    }
    colnames(exprs(x)) <- NULL
    assign(f, x, envir=env, inherits=FALSE)
  }

  return(new("cytoSet",
    frames    = env,
    phenoData = new("phenoData",
       pData = data.frame(framename=I(files)),
       varLabels = list(framename="Filename of the FCS 3.0 file")),
    colnames=cn))
}
  

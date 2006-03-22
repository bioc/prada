touchFCS <- function(path=".", file){
  ret <- NULL
  if(missing(file)){
    if(length(path)!=1 || !is.character(path))
      stop("'path' must be character of length 1")
  file <- dir(path, full.names=TRUE)
  fnames <- dir(path)
  }else{
    if(length(file)!=1 || !is.character(file))
      stop("'file' must be character of length 1")
    fnames <-  gsub(".*[/\\]", "", file)
  }
  if(length(file)>0){
    for(i in 1:length(file)){
      head <- readChar(file[i], 6)
      if(head %in% c("FCS2.0", "FCS3.0"))
        ret <- c(ret, fnames[i])
    }
  }
  return(ret)
}

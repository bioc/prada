assert.data.frame = function(x, onerow=TRUE) {
  
  if (!is.data.frame(x)) {
    if (is.character(x)) {
      stop(x)
    } else {
      stop("x is not a data.frame")
    }
  }
  
  if(onerow==TRUE && nrow(x) != 1)
    stop(paste("Expected one row in table 'x' but found", nrow(x), "\n"))

  return(TRUE)
}

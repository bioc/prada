as.all = function(x, what) {
  na = is.na(x)
  x  = do.call(paste("as.", what, sep=""), list(x))
  if(any(is.na(x)!=na))
    stop("Non-numeric values were produced in as.all")
  return(x)
}


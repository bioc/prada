## apply a statistic to the data from each well in a plate
analysePlate <-  function(x, wellcol="well", statfun, platename, plotdir=".", ...) {

  jw <- which(colnames(dat)==wellcol)
  if(length(jw)!=1)
    stop(paste("'x' must contain exactly one column with name '", wellcol,
    "'", sep=""))

  stopifnot(is.character(statfun), is.character(plotdir), is.character(platename))
  stopifnot(length(statfun)==1, length(plotdir)==1, length(platename)==1)

  y  <- NULL
  for (w in sort(unique(x[, jw]))) {
    rv <- do.call(statfun, list(x = x[x[, jw]==w, ],
                  plotdir  = outdir, ...))
    stopifnot(all(sapply(rv, length) == 1))
    crv <- sapply(rv, class)
    stopifnot(all(crv %in% c("numeric", "integer", "logical", "character")))
    
    if(getPradaPar("debug"))
      cat(".") 
    
    ## to prevent as.data.frame.character being called and converting to factor:
    for (i in which(crv=="character"))
      class(rv[[i]]) <- "vector"
    y  <- rbind(y, data.frame(platename, w, rv))
  }
  
  if(getPradaPar("debug"))
    cat("\n")
  if(is.null(names(rv)) && length(rv)==1)
    names(rv) <- statfun
  stopifnot(!is.null(names(rv)))
  
  colnames(y) <- c("platename", wellcol, names(rv))
  return(y)
}


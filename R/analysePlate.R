## apply a statistic to the data from each well in a given (expId, expRepeat)
analysePlate <-  function(eid, er, stat, pars, outdir, nrwell=96) {
  stopifnot(is.character(eid) && is.numeric(er))
  stopifnot(length(eid)==1 && length(er)==1)
  
  ## to prevent as.data.frame.character being called and converting to factor:
  platename <- paste(eid, er, sep="_")
  if(pars$debug)
    cat("analysePlate ", platename, ":", sep="")
  
  datwh <- dat[(dat$expId==eid) & (dat$expRepeat==er), ]

  ct <- estimateCrosstalk(datwh, device="png", 
        plotfile = file.path(outdir, paste(platename, "--ct.png", sep="")),
        pars     = pars)
  if(pars$debug)
    cat(" ct:", paste(names(ct), signif(ct, 2), sep=":", collapse="\t"), "\n")
  
  y  <- NULL
  for (w in 1:nrwell) {
    rv <- do.call(stat,
      list(    x      = datwh[datwh$well==w, ],
                pars  = append(pars, list(ct=ct)),
              outdir  = outdir,
             plotfile = paste(platename, "--", w, ".png", sep="")))
    crv <- sapply(rv, class)

    if(pars$debug) {
      cat(".") 
      stopifnot(all(sapply(rv, length) == 1))
      stopifnot(all(crv %in% c("numeric", "integer", "logical", "character")))
    }
    
    ## to prevent as.data.frame.character being called and converting to factor:
    for (i in which(crv=="character"))
      class(rv[[i]]) <- "vector"
    y  <- rbind(y, data.frame(I(eid), er, w, rv))
  }
  
  if(pars$debug) cat("\n")
  if(is.null(names(rv)) && length(rv)==1) names(rv) = stat
  stopifnot(!is.null(names(rv)))
  
  colnames(y) <- c("expId", "expRepeat", "well", names(rv))
  return(y)
}


estimateCrosstalk <- function(x, pars, device, plotfile=NULL)
{
  if(pars$debug)
    stopifnot(is.data.frame(x),
      "minNrCells" %in% names(pars), 
      all(c("cloneId", "well", "Field", "trsf", "dapi") %in% colnames(x)))
   
  dye       <- getDye(as.character(x$cloneId))
  thedyes   <- unique(dye)
  rv        <- numeric(length(thedyes))
  names(rv) <- thedyes

  ## background correction
  ffield <- factor(paste(x$well, x$Field, sep=":"))
  meds   <- tapply(x$trsf, ffield, median)
  trsfbg <- x$trsf - meds[as.character(ffield)]

  for (d in thedyes) {
    sel <- which(dye==d)
    if(length(sel) >= pars$minNrCells) {
      y <- data.frame(dapi=x$dapi[sel], trsfbg=trsfbg[sel])
      a <- rlm(trsfbg ~ dapi, data=y)
      rv[d] <- coef(a)["dapi"]
      
      if(!missing(device)) {
        do.call(device, args=list(file=plotfile))
        scale <- a$s
        plot(y$dapi, y$trsfbg, ylim=scale*c(-3,8), pch=".",
             main=paste(eid, "_", er, "(", dye, ") ", nrow(x), " rows",
                        " slope=", signif(coef(a)["dapi"],2),
                        " scale=", signif(scale,2), sep=""))
        for(i in (-2:2))
          abline(a=coef(a)["(Intercept)"]+scale*i,
                 b=coef(a)["dapi"],
                 col=c("red","black","grey")[abs(i)+1])
        if(device %in% c("pdf", "png", "jpeg"))
          dev.off()
      } ## if
    } ## if
  } ## for
  return(rv)
}


estimateCrosstalkWell <- function(x, do.plot=TRUE) {
  slope  <- as.numeric(NA)
  ffield <- factor(x$Field)
  nlev   <- nlevels(ffield)
  if ((nlev>=2) && (nrow(x) >= 50+nlev)) {
    a         <- rlm(x$trsf ~ x$dapi + ffield)
    stopifnot(identical(names(coef(a)),
         c("(Intercept)", "x$dapi", paste("ffield", levels(ffield)[2:nlev], sep=""))))
    slope <- coef(a)[2]
    bg    <- coef(a)[1] + c(0, coef(a)[-(1:2)])
    names(bg) <- levels(ffield)
    bg    <- bg[as.character(ffield)]  ## a value for each cell

    if(do.plot) {
      scale <- a$s
      plot(x$dapi, x$trsf-bg, ylim=scale*c(-3,8), pch=".", 
           main=paste("slope=", signif(slope, 2)))
      points(x$dapi, a$fitted.values-bg, pch=16, col="blue")
      abline(a=0, b=slope, col="red")
    }
  }
  return(slope)
}

estimateCrosstalkPlate <- function(x, pars, plotfileprefix=NULL)
{
  if(pars$debug)
    stopifnot(is.data.frame(x),
      "minNrCells" %in% names(pars), 
      all(c("cloneId", "well", "Field", "trsf", "dapi") %in% colnames(x)))
   
  dye       <- getDye(as.character(x$cloneId))
  thedyes   <- unique(dye)

  slopes <- vector(mode="list", length=length(thedyes))
  names(slopes) <- thedyes
  
  for (d in thedyes) {
    sel <- which(dye==d)
    if(length(sel) >= pars$minNrCells) {
      if(!is.null(plotfileprefix)) {
        png(width=1024, height=768, file=paste(plotfileprefix, "_", d, ".png", sep=""))
        par(mfrow=c(6, 8), mai=c(0.02,0.02,0.25,0.02))
      }
      
      y <- x[sel, ]
      slopes[[d]] <- by(y, factor(y$well), estimateCrosstalkWell,
                        do.plot=!is.null(plotfileprefix))
      
      if(!is.null(plotfileprefix))
        dev.off()
    }
  }
  rv <- sapply(slopes, shorth, na.rm=TRUE)
  
  if(!is.null(plotfileprefix)) {
    png(width=768, height=384, file=paste(plotfileprefix, "_hist.png", sep=""))
    par(mfrow=c(1,2))
    cols        <- c("#984ea3", "#4daf4a")
    names(cols) <- c("cfp",     "yfp")
    for(i in 1:length(slopes))
      hist(slopes[[i]], xlab="slopes", main=paste(thedyes[i], signif(rv[i], 2)),
           col=cols[thedyes[i]])
    dev.off()
  }
  
  return(rv)
}
     
##    scale <- a$s
##    plot(y$dapi, y$trsfbg, ylim=scale*c(-3,8), pch=".",
##        main=paste(eid, "_", er, "(", d, ") ", nrow(x), " rows",
##                   " slope=", signif(coef(a)["dapi"],2),
##                   " scale=", signif(scale,2), sep=""))
##    for(i in (-2:2))
##      abline(a=coef(a)["(Intercept)"]+scale*i,
##             b=coef(a)["dapi"],
##             col=c("red","black","grey")[abs(i)+1])
## } ## if
 

estimateCrosstalkWell <- function(x, doPlot=TRUE) {
  slope  <- as.numeric(NA)
  ffield <- factor(x$Field)
  nlev   <- nlevels(ffield)
  dkappa <- as.double(NA)
  if ((nlev>=2) && (nrow(x) >= 50+nlev)) {
    a         <- rlm(x$trsf ~ x$dapi + ffield)
    stopifnot(identical(names(coef(a)),
         c("(Intercept)", "x$dapi", paste("ffield", levels(ffield)[2:nlev], sep=""))))
    dkappa <- coef(a)[2]

    if(doPlot) {
      bg    <- coef(a)[1] + c(0, coef(a)[-(1:2)])
      names(bg) <- levels(ffield)
      bg    <- bg[as.character(ffield)]  ## a value for each cell
      
      plot(x$dapi, x$trsf-bg, ylim= a$s * c(-3,8), pch=".", 
           main=paste("kappa=", signif(dkappa, 2)))
      points(x$dapi, a$fitted.values-bg, pch=16, col="blue")
      abline(a=0, b=dkappa, col="red")
    }
  }
  return(dkappa)
}

estimateCrosstalkPlate <- function(x, plotfileprefix=NULL)
{
  stopifnot(is.data.frame(x),
            all(c("cloneId", "well", "Field", "trsf", "dapi") %in% colnames(x)))
   
  dye       <- getDye(as.character(x$cloneId))
  thedyes   <- unique(dye)

  kappas <- vector(mode="list", length=length(thedyes))
  names(kappas) <- thedyes
  
  for (d in thedyes) {
    sel <- which(dye==d)
    if(length(sel) >= getPradaPar("minNrCells")) {
      if(!is.null(plotfileprefix)) {
        png(width=1024, height=768, file=paste(plotfileprefix, "_", d, ".png", sep=""))
        par(mfrow=c(6, 8), mai=c(0.25,0.02,0.25,0.02))
      }
      
      y <- x[sel, ]
      kappas[[d]] <- by(y, factor(y$well), estimateCrosstalkWell,
                        doPlot=!is.null(plotfileprefix))
      
      if(!is.null(plotfileprefix))
        dev.off()
    }
  }
  shorths <- sapply(kappas, shorth, na.rm=TRUE)
  mads    <- sapply(kappas, mad,    na.rm=TRUE)
  
  if(!is.null(plotfileprefix)) {
    png(width=768, height=384, file=paste(plotfileprefix, "_hist.png", sep=""))
    par(mfrow=c(1,2))
    cols        <- c("#984ea3", "#4daf4a")
    names(cols) <- c("cfp",     "yfp")
    for(i in 1:length(kappas))
      hist(kappas[[i]], xlab="kappa",
           main=paste(thedyes[i], signif(shorths[i], 2), "+/-", signif(mads[i], 2)),
           col=cols[thedyes[i]])
    dev.off()
  }
  
  return(shorths)
}
     

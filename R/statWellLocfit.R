statWellLocfit = function(x, plotfileprefix, crosstalk, span) {
  stopifnot(all(c("brdu", "trsf", "dapi", "Field", "cloneId") %in% colnames(x)),
            is.numeric(crosstalk), length(crosstalk)==1, !is.na(crosstalk),
            is.numeric(span),      length(span)==1,      !is.na(span))

  ## default values
  nrcells <- nrow(x)
  trsfeff <- as.numeric(NA)
  delta   <- se.delta <- zscore <- as.numeric(NA)
  cloneId <- as.character(NA)
  alm     <- NULL
  bg <- numeric(nrow(x))
  fg <- function(x) numeric(length(x))
  
  if (nrcells>0) {
    cloneId   <- as.character(unique(x$cloneId))
    expId     <- as.character(unique(x$expId))
    expRepeat <- as.character(unique(x$expRepeat))
    expWell   <- unique(x$well)
    dye       <- as.character(unique(getDye(as.character(x$cloneId))))
    rgtau     <- getPradaPar("minRgTau")[dye]
    
    stopifnot(length(cloneId)==1, length(expId)==1, length(expRepeat)==1,
              length(expWell)==1, length(dye)==1,   !is.na(rgtau))

    mx <- data.frame(tau    = x$trsf - crosstalk * x$brdu,
                     brdu   = x$brdu,
                     ffield = factor(x$Field))

    ## check whether range of tau is large enough
    if(nrcells>40) {
      stau    <- sort(mx$tau)
      trsfeff <- (stau[nrcells-19] - stau[20]) / rgtau
    } else {
      trsfeff <- 0
    }
    if ((nlevels(mx$ffield)>=2) &&
        (nrow(mx) >= 3 + nlevels(mx$ffield)) &&
        (trsfeff >= 1) &&
        (nrcells >= getPradaPar("minNrCells"))) {
      alm   <- rlm(brdu ~  tau + I(tau*tau) + ffield, data = mx)
      sma   <- summary(alm)
      coefa <- coef(alm)[c("tau", "I(tau * tau)")]

      bg <- predict(alm, newdata=data.frame(tau=numeric(nrow(mx)), ffield=mx$ffield))
      fg <- function(x) { (coefa[1] + coefa[2] * x) * x }
      
      stopifnot(almostEqual(fg(mx$tau)+bg, fitted.values(alm)))
      stopifnot(identical(coefa, coef(sma)[c("tau", "I(tau * tau)"), "Value"]))

      ## for visualization
      lf  <- locfit.robust(x=mx$tau, y=mx$brdu, base=bg, deg=1, alpha=span)
      
      ## for slope estimation
      dlf <- locfit.robust(x=mx$tau, y=mx$brdu, base=bg, deg=1, alpha=span, deriv=1)

      tauzero  <- shorth(mx$tau)
      slp      <- preplot(dlf, newdata=tauzero, band="local")
      delta    <- slp$fit     ## point estimate
      se.delta <- slp$se.fit  ## estimated local standard deviation
      zscore   <- delta/se.delta
      stopifnot(almostEqual(tauzero, slp$xev$xev)) ## evaluation point

      do.lines <- FALSE
      if(do.lines) {
        crv <- preplot(lf, newdata=tauzero, band="local")
      }
    }

    main <- cloneId
    cat(main, paste(c("delta", "se.delta", "zscore"), signif(c(delta, se.delta, zscore), 2),
                    sep="=", collapse="\t"), "\n")

    myplot <- function(xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,...) {
      colpal <- c("#4db8a4", "#377db8", "#e31a1c")
      cols <- colorramp(colpal)(nrow(mx))[rank(mx$tau)]
      px   <- mx$tau
      py   <- mx$brdu - bg
      if(is.null(xmin)) xmin<-quantile(px, 0.02)
      if(is.null(xmax)) xmax<-quantile(px, 0.98)
      if(is.null(ymin)) ymin<-quantile(py, 0.02)
      if(is.null(ymax)) ymax<-quantile(py, 0.98)
      plot(px, py, pch=16, col=cols, xlim=c(xmin,xmax), ylim=c(ymin,ymax), main=main, 
           xlab="Signal intensity (transfection)", ylab="BrdU intensity",
           cex.lab=1.4, cex.main=1.4, ...)
      abline(v=tauzero, lwd=3, col="#808080")
      par(new=TRUE)
      plot(lf, lwd=3, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="", ylab="")
      
      if(do.lines) {
        myline <- function(slop) {
          abline(a=crv$fit-slop*tauzero, b=slop, lty=3)
        }
        myline(delta); myline(delta-se.delta); myline(delta+se.delta)
      }
    }

    myplot()
    if(!missing(plotfileprefix)) {
      fn <- paste(paste(plotfileprefix, expId, expRepeat, expWell, sep="_"))
      plotfile1 <- savetiff(paste(fn, "full", sep="_"), density=150)
    }
    myplot(xmax=tauzero+rgtau*2)
    if(!missing(plotfileprefix))
      plotfile2 <- savetiff(paste(fn, "zoom", sep="_"), density=150)
  } ## if (nrcells>0)
  
  rv = list(nrcells=nrcells, trsfeff=trsfeff,
            delta=delta*rgtau, se.delta=se.delta*rgtau, zscore=zscore,
            plotfile1=plotfile1, plotfile2=plotfile2, 
            cloneId=cloneId, plotfile=NA)
  return(rv)
}

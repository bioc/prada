statWellLocfit = function(x, plotfileprefix, crosstalk, xmax, span) {
  stopifnot(all(c("brdu", "trsf", "dapi", "Field", "cloneId") %in% colnames(x)),
            is.numeric(crosstalk), length(crosstalk)==1, !is.na(crosstalk),
            is.numeric(xmax),      length(xmax)==1,      !is.na(xmax),
            is.numeric(span),      length(span)==1,      !is.na(span))

  ## default values
  nrcells <- nrow(x)
  trsfeff <- pval <- as.numeric(NA)
  delta   <- as.numeric(NA)
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
    rgtau     <- getPradaPar("minRgTau")
    
    stopifnot(length(cloneId)==1, length(expId)==1, length(expRepeat)==1,
              length(expWell)==1, length(dye)==1,
              dye %in% names(rgtau))

    mx <- data.frame(tau    = x$trsf - crosstalk * x$brdu,
                     brdu   = x$brdu,
                     ffield = factor(x$Field))

    ## check whether range of tau is large enough
    if(nrcells>40) {
      stau    <- sort(mx$tau)
      trsfeff <- (stau[nrcells-19] - stau[20]) / rgtau[dye]
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
      
      stopifnot(almostIdentical(fg(mx$tau)+bg, fitted.values(alm)))
      stopifnot(identical(coefa, coef(sma)[c("tau", "I(tau * tau)"), "Value"]))

      ## tval  <- coef(sma)["tau", "t value"]
      ## pval1 <- pt(tval, df=sma$df[1], lower.tail=TRUE)
      ## pval2 <- pt(tval, df=sma$df[1], lower.tail=FALSE)
      ## pval  <- 2*min(pval1, pval2)

      lf  <- locfit.raw(x=mx$tau, y=mx$brdu, base=bg, deg=1, alpha=span)
      dlf <- locfit.raw(x=mx$tau, y=mx$brdu, base=bg, deg=1, alpha=span, deriv=1)
      fg  <- function(x) {predict(lf, newdata=x)}
      tauzero <- shorth(mx$tau)
      delta   <- predict(dlf, newdata=tauzero)
      cat(main, "\tdelta: ", signif(delta, 2), "\n")
    }

    
    myplot <- function(xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,...) {
      colpal <- c("#4db8a4", "#377db8", "#e31a1c")
      main <- paste(paste(expId, expRepeat, expWell, sep="_"), cloneId)
      cols <- colorramp(colpal)(nrow(mx))[rank(mx$tau)]
      px   <- mx$tau
      py   <- mx$brdu - bg
      if(is.null(xmin)) xmin<-quantile(px, 0.02)
      if(is.null(xmax)) xmax<-quantile(px, 0.98)
      if(is.null(ymin)) ymin<-quantile(py, 0.02)
      if(is.null(ymax)) ymax<-quantile(py, 0.98)
      plot(px, py, pch=16, col=cols, main=main, xlim=c(xmin,xmax), ylim=c(ymin,ymax),
           xlab="Signal intensity (transfection)", ylab="BrdU intensity",
           cex.lab=1.4, cex.main=1.4,...)
      
      px <- seq(from=min(mx$tau), to=max(mx$tau), length=100)
      lines(px, fg(px), lwd=3, col="black")
      abline(v=tauzero, lwd=3, col="#808080")
      
    }

    dev.set(dev.next())
    myplot()
    if(!missing(plotfileprefix)) {
      fn <- paste(paste(plotfileprefix, expId, expRepeat, expWell, sep="_"))
      savetiff(paste(fn, "full", sep="_"), density=150)
    }
    
    dev.set(dev.next())
    myplot(xmax=xmax)
    
    if(!missing(plotfileprefix))
      savetiff(paste(fn, "zoom", sep="_"), density=150)
  } ## if (nrcells>0)
  
  rv = list(nrcells=nrcells, trsfeff=trsfeff, delta=delta, 
             pval=pval, cloneId=cloneId, plotfile=NA)
  return(rv)
}

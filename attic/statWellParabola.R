statWellParabola = function(x, plotfile=as.character(NA), crosstalk) {
  stopifnot(all(c("brdu", "trsf", "dapi", "Field", "cloneId") %in% colnames(x)))

  ## default values
  nrcells <- nrow(x)
  trsfeff <- pval <- as.numeric(NA)
  delta   <- c(as.numeric(NA), as.numeric(NA))
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
              dye %in% names(crosstalk),
              dye %in% names(rgtau))

    xtau  <- x$trsf - crosstalk[dye] * x$brdu
    xtau0 <- shorth(xtau)  ## 0, to make coefficients more interpretable?
    
    mx <- data.frame(tau    = xtau-xtau0,
                     brdu   = x$brdu,
                     ffield = factor(x$Field))

    ## check whether range of tau is large enough
    if(nrcells>40) {
      stau    <- sort(xtau)
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

      delta.raw <- coef(alm)[c("tau", "I(tau * tau)")]
      delta     <- delta.raw * rgtau[dye]
      
      stopifnot(identical(delta.raw, coef(sma)[c("tau", "I(tau * tau)"), "Value"]))

      tval  <- coef(sma)["tau", "t value"]
      pval1 <- pt(tval, df=sma$df[1], lower.tail=TRUE)
      pval2 <- pt(tval, df=sma$df[1], lower.tail=FALSE)
      pval  <- 2*min(pval1, pval2)

      bg <- predict(alm, newdata=data.frame(tau=numeric(nrow(mx)), ffield=mx$ffield))
      fg <- function(x) { (delta.raw[1] + delta.raw[2] * (x-xtau0)) * (x-xtau0) }
      stopifnot(almostEqual(fg(xtau)+bg, fitted.values(alm)))
    }

    if(!is.na(plotfile)) {
      png(file = plotfile, width=768, height=512)
      layout(matrix(c(1,3,2,4), nrow=2, ncol=2), widths=c(1,1), heights=c(3,2))
      colpal <- c("#4db8a4", "#377db8", "#e31a1c")

      myplot <- function(...) {
        cols <- colorramp(colpal)(nrow(mx))[rank(mx$tau)]
        plot(xtau, mx$brdu - bg, pch=16, col=cols,
             xlab="transfection", ylab="BrdU", ...)

        px <- seq(from=min(xtau), to=max(xtau), length=100)
        lines(px, fg(px), lwd=2, col="black")
        abline(v=xtau0, col="#808080")
      }
      myplot(main=paste(paste(expId, expRepeat, expWell, sep="_"), cloneId),
             xlim=min(xtau)+c(0,3)*rgtau[dye])
      myplot(main=paste("delta=", paste(signif(delta,2), collapse=" "),
                        " p=", format.pval(pval,2), sep=""))

      if(!is.null(alm)) {
        mai <- par("mai"); mai[3] <- 0; par(mai=mai)
        plot(rank(xtau), residuals(alm), pch=16, xlab="rank(transfection)")
        abline(h=0, col="#e31a1c")
        plot(bg, pch=".", col="#4db8a4", xlab="cells", ylab="BrdU bg")
      }
      dev.off()
    } ## if (!missing(plotfile))
  } ## if (nrcells>0)
  
  rv = list(nrcells=nrcells, trsfeff=trsfeff, delta1=delta[1], delta2=delta[2],
             pval=pval, cloneId=cloneId, plotfile=plotfile)
  return(rv)
}

statWellLocfit = function(x, span, plotwhat="nothing", plotdir=".", plotfile, ...) {
  
  stopifnot(all(c("brdu", "trsf", "dapi", "Field", "cloneId") %in% colnames(x)),
            is.numeric(span),      length(span)==1,      !is.na(span),
            is.character(plotwhat), length(plotwhat)==1,
            is.character(plotdir),  length(plotdir)==1)
  
  nrcells <- nrow(x)
  if(nrcells==0) {
    num <- as.numeric(NA)
    chr <- as.character(NA)
    return(list(cloneId=chr, nrcells=num, trsfeff=num,
                delta=num, se.delta=num, zscore=num, plotfile=chr))
  }
  
  cloneId   <- as.character(unique(x$cloneId))
  expId     <- as.character(unique(x$expId))
  expRepeat <- as.character(unique(x$expRepeat))
  expWell   <- unique(x$well)
  dye       <- unique(getDye(as.character(x$cloneId)))
  stopifnot(length(cloneId)==1, length(expId)==1, length(expRepeat)==1,
            length(expWell)==1, length(dye)==1,
            !is.na(cloneId), !is.na(expId), !is.na(expRepeat),
            !is.na(expWell), !is.na(dye))

  rgtau   <- getPradaPar("minRgTau")[dye]
  stopifnot(!is.na(rgtau))

  tau     <- x$trsf
  tauzero <- shorth(tau)
  ffield  <- factor(x$Field)

  ## default values:
  delta <- se.delta <- zscore <- trsfeff <- as.numeric(NA)
  ybg   <- x$brdu
  lcft <- resi <- NULL
  
  ## check whether range of tau ("transfection efficiency") is large enough
  if(nrcells>40) {
    stau    <- sort(tau)
    nrctb   <- getPradaPar("nrCellsTopBottom")
    trsfeff <- (stau[nrcells-nrctb+1] - stau[nrctb]) / rgtau
  } else {
    trsfeff <- 0
  }
  if ((nlevels(ffield)>=2) &&
      (nrcells >= 3 + nlevels(ffield)) &&
      (trsfeff >= 1) &&
      (nrcells >= getPradaPar("minNrCells"))) {

    bg     <- numeric(nrcells)

    if(TRUE) {
      ## with background
      sc.dbg <- sc.res.old <- Inf
      niter  <- 0
      tol    <- 1
      while((sc.dbg > tol) && (niter < 10)) {
        lcft     <- locfit.robust(x=tau, y=x$brdu - bg, deg=1, alpha=span, maxk=512)
        resi     <- residuals(lcft)
        sc.res   <- mad(resi)
        ## stopifnot(sc.res <= sc.res.old)  ## should monotonously decrease ?
        sc.res.old <- sc.res
        bgff     <- tapply(resi, ffield, median)
        dbg      <- bgff[as.character(ffield)]
        sc.dbg   <- diff(range(bgff))
        bg       <- bg + dbg
        bg       <- bg - mean(bg)
        stopifnot(!any(is.na(bg)))
        niter    <- niter+1
      }
    } else {
      ## without background
      lcft <- locfit.robust(x=tau, y=x$brdu, deg=1, alpha=span, maxk=512)
      resi <- residuals(lcft)
    }
    
    ybg <- x$brdu - bg 
    ## for slope estimation
    dlcft <- locfit.robust(x=tau, y=ybg, deg=1, alpha=span, deriv=1, maxk=512)

    slp      <- preplot(dlcft, newdata=tauzero, band="local")
    delta    <- slp$fit*rgtau     ## point estimate
    se.delta <- slp$se.fit*rgtau  ## estimated local standard deviation
    zscore   <- delta/se.delta
    stopifnot(all(abs(tauzero-slp$xev$xev)< 1e-6)) ## at evaluation point
    
  } ## if (enough transfections efficiency, nrcells, etc.)
  
  myplot <- function(qxmin=0, qxmax=1, qymin=0, qymax=1, lxmax, ...) {
    colpal <- c("#4db8a4", "#377db8", "#e31a1c")
    cols <- colorRampPalette(colpal)(nrcells)[rank(tau)]
    px   <- tau
    py   <- ybg
    xlim <- quantile(px, c(qxmin, qxmax))
    ylim <- quantile(py, c(qymin, qymax))
    plot(px, py, pch=20, col=cols,
         xlim=xlim, ylim=ylim, xlab=getPradaPar("xlab")[dye], ylab=getPradaPar("ylab"),
         cex.lab=1.4, cex.main=1.4, ...)
    abline(v=tauzero, lwd=3, col="#808080")
    if(!is.null(lcft)) {
      if(missing(lxmax))
        lxmax <- sort(px)[max(1, length(px)-10)]
      px <- seq(xlim[1], lxmax, len=50)
      lines(px, predict(lcft, newdata=px), lwd=3)
    }
  }

  if (missing(plotfile))
    plotfile <- paste(expId, expRepeat, expWell, sep="_")

  switch(plotwhat,
     nothing = {},
     screen  = {
       myplot(...) 
       plotfile <- as.character(NA)
     },
     figscp  = {
       plotfile <- paste(plotfile, ".pdf", sep="")
       pdf(file=file.path(plotdir, plotfile), width=7, height=7)
       myplot(...)
       dev.off()
       cat(cloneId, paste(c("delta", "se.delta", "zscore", "niter", "sc.dbg"),
                          signif(c(delta, se.delta, zscore, niter, sc.dbg), 2),
                          sep="=", collapse="\t"), "\n")
     },
     mkpp    = {
       main1 <- paste(plotfile, paste(mapId(cloneId), "(low range)"))
       main2 <- paste("delta=", signif(delta, 2), "  se.delta=", signif(se.delta, 2),
                      "  z=", signif(zscore, 2), "  (full range)", sep="")
       plotfile <- paste(plotfile, ".png", sep="")
       png(file=file.path(plotdir, plotfile), width=1024, height=768)
       layout(matrix(c(1,3,2,4), nrow=2, ncol=2), widths=c(1,1), heights=c(2,1))
       myplot(main=main1, qxmax=0.95, ...)
       myplot(main=main2, ...)
       if(!is.null(resi)) {
         plot(rank(tau), resi, pch=16, xlab="rank(x)", ylab="residuals")
         abline(h=0, col="red")
         plot(bg, pch=".", col="#4db8a4", xlab="cells", ylab="BrdU bg")
       }
       dev.off()
     },
     stop(paste("Unknown value of 'plotwhat':", plotwhat))
  )
  names(plotfile) <- "plotfile"
  rv <- append(list(
    cloneId=cloneId, nrcells=nrcells, trsfeff=trsfeff,
    delta=delta, se.delta=se.delta, zscore=zscore),
    plotfile)
  
  return(rv)
}

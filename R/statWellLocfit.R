statWellLocfit = function(x, plotwhat="nothing", plotdir=".", crosstalk, span) {
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
  dye       <- as.character(unique(getDye(as.character(x$cloneId))))
  rgtau     <- getPradaPar("minRgTau")[dye]

  if(length(crosstalk)>1)
    crosstalk <- crosstalk[dye]
  stopifnot(is.numeric(crosstalk), !is.na(crosstalk))
  
  stopifnot(length(cloneId)==1, length(expId)==1, length(expRepeat)==1,
            length(expWell)==1, length(dye)==1,   !is.na(rgtau))

  tau     <- x$trsf - crosstalk * x$brdu
  tauzero <- shorth(tau)
  ffield  <- factor(x$Field)

  ## default values:
  delta <- se.delta <- zscore <- trsfeff <- as.numeric(NA)
  ybg   <- x$brdu
  lcft <- resi <- NULL
  
  ## check whether range of tau ("transfection efficiency") is large enough
  if(nrcells>40) {
    stau    <- sort(tau)
    trsfeff <- (stau[nrcells-19] - stau[20]) / rgtau
  } else {
    trsfeff <- 0
  }
  if ((nlevels(ffield)>=2) &&
      (nrcells >= 3 + nlevels(ffield)) &&
      (trsfeff >= 1) &&
      (nrcells >= getPradaPar("minNrCells"))) {

    ## for location estimation (including background)
    bg     <- numeric(nrcells)
    sc.dbg <- sc.res.old <- Inf
    niter  <- 0
    tol    <- 1
    while(sc.dbg > tol) {
      lcft     <- locfit.robust(x=tau, y=x$brdu - bg, deg=1, alpha=span, maxk=512)
      resi     <- residuals(lcft)
      sc.res   <- mad(resi)
      ## stopifnot(sc.res <= sc.res.old)  ## should monotonously decrease ?
      sc.res.old <- sc.res
      bgff       <- tapply(resi, ffield, median)
      dbg      <- bgff[as.character(ffield)]
      sc.dbg   <- diff(range(bgff))
      bg       <- bg + dbg
      bg       <- bg - median(bg)
      stopifnot(!any(is.na(bg)))
      niter    <- niter+1
    }
    ybg <- x$brdu - bg 

    ## for slope estimation
    dlcft <- locfit.robust(x=tau, y=ybg, deg=1, alpha=span, deriv=1, maxk=512)

    slp      <- preplot(dlcft, newdata=tauzero, band="local")
    delta    <- slp$fit     ## point estimate
    se.delta <- slp$se.fit  ## estimated local standard deviation
    zscore   <- delta/se.delta
    stopifnot(almostEqual(tauzero, slp$xev$xev)) ## evaluation point
    
  } ## if (enough transfections efficiency, nrcells, etc.)
  
  myplot <- function(xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL,...) {
    colpal <- c("#4db8a4", "#377db8", "#e31a1c")
    cols <- colorramp(colpal)(nrcells)[rank(tau)]
    px   <- tau
    py   <- ybg
    if(is.null(xmin)) xmin <- quantile(px, 0.01)
    if(is.null(xmax)) xmax <- quantile(px, 0.99)
    if(is.null(ymin)) ymin <- quantile(py, 0.01)
    if(is.null(ymax)) ymax <- quantile(py, 0.99)
    plot(px, py, pch=16, col=cols, xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
         xlab="Signal intensity (transfection)", ylab="BrdU intensity",
         cex.lab=1.4, cex.main=1.4, ...)
    abline(v=tauzero, lwd=3, col="#808080")

    if(!is.null(lcft)) {
      px <- seq(xmin, xmax, len=50)
      lines(px, predict(lcft, newdata=px), lwd=3)
    }
  }

  pfn <- paste(expId, expRepeat, expWell, sep="_")
  
  switch(plotwhat,
     nothing = {
       plotfile <- as.character(NA)
     },
     figscp  = {
       par(mfrow=c(1,1))
       myplot(main=cloneId)
       f1 <- savetiff(paste(pfn, "full", sep="_"), density=150, dir=plotdir)
       myplot(main=cloneId, xmax=tauzero+rgtau*2)
       f2 <- savetiff(paste(pfn, "zoom", sep="_"), density=150, dir=plotdir)
       plotfile <- c(f1, f2)
       names(plotfile) <- c("plotfull", "plotzoom")
       
       cat(cloneId, paste(c("delta", "se.delta", "zscore", "niter", "sc.dbg"),
                          signif(c(delta, se.delta, zscore, niter, sc.dbg), 2),
                          sep="=", collapse="\t"), "\n")
     },
     mkpp    = {
       plotfile <- paste(pfn, ".png", sep="")
       names(plotfile) <- "plotfile"
       png(file=file.path(plotdir, plotfile), width=1024, height=768)
       layout(matrix(c(1,3,2,4), nrow=2, ncol=2), widths=c(1,1), heights=c(2,1))
       myplot(main=paste(pfn, cloneId))
       myplot(main=paste("delta=", signif(delta, 2), "+/-", signif(2*se.delta, 2), " z=", signif(zscore, 2), sep=""),
              xmax=tauzero+rgtau*2)
       if(!is.null(resi)) {
         plot(rank(tau), resi, pch=16, xlab="rank(x)", ylab="residuals")
         abline(h=0, col="red")
         plot(bg, pch=".", col="#4db8a4", xlab="cells", ylab="BrdU bg")
       }
       dev.off()
     },
     stop(paste("Unknown value of 'plotwhat':", plotwhat))
  )
  
  rv <- append(list(
    cloneId=cloneId, nrcells=nrcells, trsfeff=trsfeff,
    delta=delta*rgtau, se.delta=se.delta*rgtau, zscore=zscore),
    plotfile)
  
  return(rv)
}

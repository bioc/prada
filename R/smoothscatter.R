smoothScatter <- function(x, y=NA, nrpoints=100,
                           colrange=c("white", "#231d1d"),
                           nbin=200, bandwidth, ...){
  
    # check correct input
    if (!is.matrix(x) || ncol(x)<2)
      x <- cbind(x, y)
    if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) > 1) )
           stop("'nrpoints' should be integer or inf")
    nrpoints <- floor(nrpoints)
    if (!is.character(colrange))
      stop("colrange must be character vector containing valid color identifiers")
    if (!is.numeric(nbin) || (!length(nbin) %in% 1:2))
      stop("'nbin' must be numeric of length 1 or 2")
    if (length(nbin) == 1)
      nbin <- c(nbin, nbin)

    if (missing(bandwidth)) {
      bandwidth <- diff(apply(x, 2, quantile, probs=c(0.05, 0.95))) / 10
      } else {
      if(!is.numeric(bandwidth))
        stop("'nbin' and 'bandwidth' must be numeric")
     }
                                        # create density map
    map  <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth)
    xm   <- map$x1
    ym   <- map$x2
    dens <- map$fhat
    m    <- nrow(x)
  
    # plot color image
    cols <- colorramp(colrange) (256)
    image(xm, ym, z=dens, col=cols, ...)
    box()

    # plot selection of dots
    if (nrpoints!=0){
      ## we assume that map$x1 and map$x2 go linearly from
      ## their first to their last value in nbin steps
      stopifnot(length(xm)==nbin[1], length(ym)==nbin[2], all(dim(dens)==nbin))
      ixm <- round((x[,1]-xm[1])/(xm[nbin[1]]-xm[1])*(nbin[1]-1))
      iym <- round((x[,2]-ym[1])/(ym[nbin[2]]-ym[1])*(nbin[2]-1))
      idens <- dens[1+iym*nbin[1]+ixm]
      nrpoints <- min(m, nrpoints)
      sel <- order(idens, decreasing=FALSE)[1:nrpoints]
      points(x[sel,1:2], pch=".", col="black")
    }
}

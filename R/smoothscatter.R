smoothScatter <- function(x, y=NA, nrpoints=100,
                           colrange=c("white", "#231d1d"),
                           nbin=200, bandwidth=20, ...){
  
    # check correct input
    if (!is.matrix(x) || ncol(x)<2)
      x <- cbind(x, y)
    if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) > 1) )
           stop("'nrpoints' should be integer or inf")
    nrpoints <- floor(nrpoints)
    if (!is.character(colrange))
      stop("colrange must be character vector containing valid color identifiers")
    if (!is.numeric(nbin) | !is.numeric(bandwidth))
      stop("'nbin' and 'bandwidth' must be numeric")
    if ( length(nbin) == 1 )
        nbin <- c(nbin, nbin)

    # create density map
    map  <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth)
    xm    <- map$x1
    ym    <- map$x2
    dens <- map$fhat
    m    <- nrow(x)
  
    # plot color image
    cols <- colorramp(colrange) (256)
    image(xm, ym, z=dens, col=cols, ...)
    box()
    
    # plot selection of dots
    if (nrpoints!=0){
      nrpoints <- min(m, nrpoints)
      sel <- sample(1:m, size=nrpoints)
      points(x[sel,1:2], pch=".", col="black")
    }
}

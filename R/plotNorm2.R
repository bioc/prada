plotNorm2 <- function(fn, colrange=c("gray82", "blue"), center=TRUE,
                      selection=FALSE, ellipse=FALSE, pch=20, cex=1, ...) {

  require(geneplotter)
  
  if(!is.list(fn) || any(names(fn)!=c("mu", "S", "p", "sel", "scalefac", "data")))
    stop("Parameter 'fn' must be a list with elements mu, S, p, sel, scalefac, data")
  for(i in c("center", "selection", "ellipse"))
    if(!is.logical(get(i)))
      stop(paste("Parameter ', i, ' must be logical.", sep=""))
  if(!is.character(colrange))
    stop("Parameter 'colrange' must be a character vector")

  nrcolors <- 256
  palette  <- colorRampPalette(colrange)(nrcolors)
  xrange   <- range(fn$p, na.rm=TRUE)
  z2icol   <- function(z) {
    res = round((z-xrange[1])/diff(xrange)*(nrcolors-1))+1
    res[res > nrcolors] <- nrcolors
    res[res < 1       ] <- 1
    return(res)
  }
  
  cols <- palette[z2icol(fn$p)]

  ## produce plot
  if(ellipse){
    plot(fn$data, type="n", pch=pch, cex=cex, ...)
    ell <- .addEllipse(S=fn$S, mu=fn$mu, rad=fn$scalefac)
    points(fn$data, col=cols, pch=pch, cex=cex)
    lines(ell)
  }
  else plot(fn$data, col=cols, pch=pch, cex=cex, ...)

  if(center)
    points(fn$mu[1], fn$mu[2], pch=4, col="red", cex=2)

  if(selection)
    points(fn$data[!fn$sel,], pch=".", col="red")
}

.addEllipse <- function (S, mu, rad, ...){

  stopifnot((is.matrix(S) && dim(S)==c(2,2)), is.numeric(mu), length(mu)==2,
            is.numeric(rad), length(rad)==1)

# add ellipse  
  eigen <- eigen(S)
  phi   <- seq(0, 2*pi, len = 180)
  phi0  <- -atan2(abs(eigen$vector[1,2]), abs(eigen$vector[1,1]))
  xc    <- rad * sqrt(eigen$value[1]) * cos(phi)
  yc    <- rad * sqrt(eigen$value[2]) * sin(phi)
  xc1 <- mu[1] + xc*cos(phi0) + yc*sin(phi0)
  yc1 <- mu[2] - xc*sin(phi0) + yc*cos(phi0)
  polygon(x=xc1, y=yc1, col="#fffab1")
  ell <- cbind(xc1, yc1)
  return(invisible(ell))
}


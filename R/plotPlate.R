## plot a statistic
## x must be a dataframe with columns "well" and "stat"
plotPlate  = function(x, stat, main, zlim=NULL, col, device, plotfile=NULL, width,
                      na.action="zero") {
  nrwell = 96 
  sizex  = 12
  sizey  = 8
  js = which(colnames(x)==stat)
  if(length(js)!=1)
    stop(paste("'x' must contain exactly one column with name '", stat, "'\n", sep=""))
  jw = which(colnames(x)=="well")
  if(length(jw)!=1)
    stop("'x' must contain exactly one column with name 'well'\n")
  if(any((x[,jw]<1) | (x[,jw]>nrwell)))
    stop(paste("'x$well' must be between 1 and", nrwell, "\n"))
  if(any(duplicated(x[,jw]) | is.na(x[,jw])))
    stop("'x$well' must not contain NA or duplicates.\n")
  if(is.null(zlim))
    zlim=range(x[,js], na.rm=TRUE)
  if(!is.numeric(zlim) || length(zlim)!=2 || any(is.na(zlim)))
    stop("'zlim' mist be a numeric vector of length 2 with no NA.")

  ## user coordinates: x=(-0.5...13.5), y=(-0.5...9.5)
  xlim   = c(-0.5, sizex+1.5)
  ylim   = c(-0.5, sizey+1.5)
  colbarwid = 0.3
  fw     = diff(xlim)+colbarwid
  fh     = diff(ylim)
  height = width/fw * fh
  do.call(device, args=list(file=plotfile, width=width, height=height))
  layout(matrix(1:2, ncol=2), widths=c(diff(xlim), colbarwid), heights=1)

  ## device coordinates
  u2px = function(x) (x-xlim[1]) / fw * width
  u2py = function(y) (y-ylim[1]) / fh * height

  par(mai=c(0,0,0,0))
  cex = 1.5
  plot(x=0, y=0, type="n", bty="n", xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlim=xlim, ylim=ylim)
  text((1:sizex), 0, paste(1:sizex), adj=c(0.5,0), cex=cex)
  text(0, (sizey:1), c("A","B","C","D","E","F","G","H"), adj=c(0, 0.5), cex=cex)
  if(!missing(main))
    text((sizex+1)/2, sizey+1, main, adj=c(0.5, 1), cex=cex)
  
  nrcolors   = 256
  thepalette = colorramp(col)(nrcolors)

  # the mapping from values to color indices
  z2icol = function(z) {
    res = round((z-zlim[1])/diff(zlim)*(nrcolors-1))+1
    res[res > nrcolors] = nrcolors
    res[res < 1       ] = 1
    return(res)
  }
  icol2z = function(i) {
    (i-1)/(nrcolors-1)*diff(zlim)+zlim[1]
  }
  stopifnot(all(z2icol(icol2z(1:nrcolors))==1:nrcolors))
  circcol = thepalette[z2icol(x[,js])]

  ## circles
  radius = 0.45
  xc = radius*cos(seq(0, 2*pi, len=73))
  yc = radius*sin(seq(0, 2*pi, len=73))
  x0 = 1     + ((x[,jw]-1)  %% sizex)
  y0 = sizey - ((x[,jw]-1) %/% sizex)

  switch(na.action,
         zero = {
           circcol[is.na(circcol)] <- thepalette[z2icol(0)]
           wh <- 1:nrow(x)
         }, 
         omit = {
           wh <- which(!is.na(circcol))
         },
         stop(paste("Invalid value of 'na.action':", na.action))
  )
  
  for (i in wh) {
    polygon(x = x0[i]+xc,
            y = y0[i]+yc,
            col=circcol[i])
  }
  
  xmin = 0.5
  xmax = sizex + 0.5
  ymin = 0.5
  ymax = sizey + 0.5
  polygon(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin))
        
  par(mai=0.5*c(1,0,1,0))
  image(x=0, y=icol2z(1:nrcolors), z=matrix(1:nrcolors, nrow=1),
        col = thepalette, xaxt="n")

  if(device %in% c("pdf", "png", "jpeg"))
    dev.off()

  x0 = 1 + (wh-1) %%  sizex
  y0 = 1 + (wh-1) %/% sizex
  dx = dy = 0.4
  x1 = u2px(x0-dx)
  x2 = u2px(x0+dx)
  y1 = u2py(y0-dy)
  y2 = u2py(y0+dy)

  res = list(which=wh, coord=floor(cbind(x1, y1, x2, y2) + 0.5))
  
  class(res) = "plotPlateResult"
  return(res)
}

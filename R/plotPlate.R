plotPlate <- function (x, nrow = 8, ncol = 12, ind = 1: (ncol*nrow), main,
                       xrange, col, device="x11", file, width, na.action = "zero",
                       desc = as.character(c(NA, NA)), char){

  require(geneplotter)
  require(grid)
  ## validate entries
  if (!is.numeric(ncol) || length(ncol) != 1)
    stop("'ncol' must be a numeric vector of length 1")
  if (!is.numeric(nrow) || length(nrow) != 1)
    stop("'nrow' must be a numeric vector of length 1")
  nrwell <- ncol * nrow
  if ((length(ind) != length(x)) | (max(ind, na.rm = TRUE) > nrwell))
    stop("'ind' must be vector of unique indices for vector 'x'")
  y <- rep(NA, nrwell)
  xorig <-  x
  y[ind] <- x
  x <- y
  if (!is.numeric(x) || !is.vector(x) || length(x) != nrwell)
    stop(paste("'x' must be a numeric vector of length 'ncol*nrow'.",
               "\nYou might want to include indices for missing wells."))
  if (missing(xrange))
    xrange = range(x, na.rm = TRUE)
  if (!is.numeric(xrange) || length(xrange) != 2 || any(is.na(xrange)))
    stop("'xrange' must be a numeric vector of length 2 with no NA.")
  if (!is.character(desc) || length(desc) != 2)
        stop("'desc' must be a character vector of length 2")
  if(!missing(char)){
    if (!is.vector(char) || length(char) != length(ind) ||
        !all(nchar(char[which(!is.na(char))])==1))
      stop(paste("\n'char' must be a numeric vector of length 'ncol*nrow'",
                 "\nwith vector items nchar==1 or 'NA'.",
                 "\nYou might want to include indices for missing wells."))
    y[ind] <- char
    char <- y
  }

  ## device & font size 
  height=((width-0.01*width)/(ncol+1))*(nrow+1)
  args <- list(width = width, height = height)
    if (!is.character(device))
      stop("'device' must be a character")
    if (!missing(file))
      args <- append(args, file)
    if(device %in% c("X11", "x11", "windows", "quartz", "gnome", "GTK"))
      if(names(dev.cur())!="null device")
        dev.off()
    do.call(device, args = args)
  font <- ceiling(12 * width/9)
  if (device %in% c("png", "jpeg"))
    font <- ceiling(12 * width/900)
  
  ## units 2 pixel 
  xlim = c(0, ncol + 1)
  ylim = c(0, nrow + 1)
  fw = diff(xlim)/0.9
  fh = diff(ylim)/0.9
  u2px = function(x) (x - xlim[1])/fw * width
  u2py = function(y) (y - ylim[1])/fh * height

  ## create false colors 
  nrcolors = 256
  thepalette = colorRampPalette(col)(nrcolors)
  z2icol <- function(z) {
    res = round((z - xrange[1])/diff(xrange) * (nrcolors -
      1)) + 1
    res[res > nrcolors] = nrcolors
    res[res < 1] = 1
    return(res)
  }
  icol2z <- function(i) {
    (i - 1)/(nrcolors - 1) * diff(xrange) + xrange[1]
  }
  stopifnot(all(z2icol(icol2z(1:nrcolors)) == 1:nrcolors))
  circcol <- thepalette[z2icol(x)]
  
  ## deal with NA values
  wh <- (1:nrwell)[ind]
  whorig <- seq(along=xorig)
  switch(na.action, zero = {
    circcol[is.na(circcol)] <- thepalette[z2icol(0)]
  }, omit = {
    wh <- which(!is.na(circcol))
    whorig <- which(!is.na(xorig))
  }, xout = {
    nawell <- which(is.na(circcol))
    sel <- nawell[which((nawell %in% ind))]
    circcol[is.na(circcol)] <- "lightgray"
  }, stop(paste("Invalid value of 'na.action':", na.action)))
  
  
  ## create grid graphic
  vp1 <- viewport(width=0.9, x=0, just="left") 
  vp2 <- viewport(width=0.1, x=0.9, just="left")
  pushViewport(vp2)
  vp3 <- viewport(height=0.85, width=0.8)
  pushViewport(vp3)
  if (any(!is.na(desc)))
    grid.text(y=c(0.95, 0.05), x=0.1, just="left", desc, gp=gpar(fontsize=font,
                    cex=1.4, fontface="bold", col=thepalette[c(length(thepalette), 1)]))
  vp4 <- viewport(height=0.8, width=0.1, yscale=c(xrange), xscale=c(0,1),
                  x=0.1, just="left")

  pushViewport(vp4)
  cr <- colorRampPalette(col) 
  cols <- thepalette[floor(seq(1, nrcolors, length=50))]
  for(i in 1:50)
    grid.rect(y=0+i/50, height=1/50,gp=gpar(fill=cols[i], col=cols[i]),
              just="top")
  grid.rect()
  at <- signif(seq(xrange[1], xrange[2], length=6)[2:5],2)
  grid.yaxis(at=at, gp=gpar(fontsize=font, cex=1), main=FALSE, label=FALSE)
  grid.text(x=unit(3.5, "native"), y=unit(at, "native"), at, rot=90,
            gp=gpar(fontsize=font, cex=1))
  popViewport(3)
  pushViewport(vp1)
  vp5 <- viewport(height=0.9, y=0, just="bottom", xscale=c(0, ncol+1),
                  yscale=c(0, nrow+1))
  pushViewport(vp5)
  grid.rect(width=unit(1-(1/(ncol+1)), "npc"), height=unit(1-(1/(nrow+1)),
            "npc"), x=unit(1/(ncol+1), "npc"), y=unit(1/(nrow+1), "npc"),
            just=c("left","bottom"))
  radius = 0.5  
  x0 = (1 + (0:(nrwell - 1))%%ncol) + radius
  y0 = (nrow - (0:(nrwell - 1))%/%ncol) + radius
  grid.circle(x=unit(x0[wh], "native"), y=unit(y0[wh], "native"),
              r=unit(radius-0.02, "native"),
              gp=gpar(fill=circcol[wh]))
  grid.text(x=unit(unique(x0), "native"), y=unit(0.9, "native"), 1:ncol, just="top",
            gp=gpar(fontsize=font, cex=1.6, fontface="bold"))
  grid.text(x=unit(0.9, "native"), y=unit(unique(y0), "native"), LETTERS[1:nrow],
            gp=gpar(fontsize=font, cex=1.6, fontface="bold"), just="right")
  if(na.action=="xout")
    for (i in sel){
      grid.lines(unit(c(x0[i]-0.39, x0[i]+0.39), "native"), unit(c(y0[i]-0.39,
                 y0[i]+0.39), "native"), gp=gpar(lwd=2, col="darkgray")) 
      grid.lines(unit(c(x0[i]-0.39, x0[i]+0.39), "native"), unit(c(y0[i]+0.39,
                 y0[i]-0.39), "native"), gp=gpar(lwd=2, col="darkgray"))
    }
  if(!missing(char))
    for (i in which(!is.na(char)))
      grid.text(x=unit(x0[i], "native"), y=unit(y0[i], "native"), char[i],
                gp=gpar(fontsize=font, cex=1.8))
  popViewport(2)
  if (!missing(main)){
    vp6 <- viewport(height=0.1, y=0.9, just="bottom")
    pushViewport(vp6)
    grid.text(main, gp=gpar(fontsize=font, cex=1.8, fontface="bold"))
    popViewport()
  }

  ## imageMap coordinates     
  if (device %in% c("pdf", "png", "jpeg"))
    dev.off()
  x0 = 1.5 + (wh - 1)%%ncol
  y0 = 0.1 * diff(ylim) + 0.6 + (wh - 1)%/%ncol
  dx = dy = 0.4
  x1 = u2px(x0 - dx)
  x2 = u2px(x0 + dx)
  y1 = u2py(y0 - dy)
  y2 = u2py(y0 + dy)

  res <- list(which = whorig, coord = floor(cbind(x1, y1, x2, y2) +
                                0.5), height = args$height, width = args$width)
  return(invisible(res))
}

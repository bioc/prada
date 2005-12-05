plotPlateGrid <- function (x, gridCall="circle", callArgs=NULL, nrow = 8, ncol = 12,
                            ind = 1: (ncol*nrow), main, xrange, col, device,
                            file, width, na.action = "zero",
                            desc = as.character(c(NA, NA)), char, legend=TRUE){

  require(geneplotter)
  require(grid)

  
  ## default plotting function ##
  circle <- function(data, col){
    grid.circle(x=0.5, y=0.5, r=0.45, gp=gpar(fill=col))
  }
  default <- ifelse(gridCall=="circle", TRUE, FALSE)
  
  ## validate parameters ##
  ## ncol
  if (!is.numeric(ncol) || length(ncol) != 1)   
    stop("'ncol' must be a numeric vector of length 1")
  ## nrow
  if (!is.numeric(nrow) || length(nrow) != 1)
    stop("'nrow' must be a numeric vector of length 1")
  nrwell <- ncol * nrow
  
  ## device
  if(!missing(device))
    if (!is.character(device) || length(device)!=1)
      stop("'device' must be a character vector of length 1")
  
  ## char
  info <- character(nrwell)
  if(!missing(char)){
    if (!is.vector(char) || length(char) != length(ind) ||
        !all(nchar(char[which(!is.na(char))])==1))
      stop(paste("\n'char' must be a  vector of length 'ncol*nrow'",
                 "\nor of length equal to inf with vector items nchar==1",
                 "or 'NA'.\nYou might want to include indices for ",
                 "missing wells."))
    char[is.na(char)] <- ""
    info[ind] <- char
  }
  ## xrange
  if(missing(xrange))
    xrange = range(x, na.rm = TRUE)
  if (!is.numeric(xrange) || length(xrange) != 2 || any(is.na(xrange)))
    stop("'xrange' must be a numeric vector of length 2 with no NA.")
  ## legend and desc
  if(!is.logical(legend) || length(legend)!=1)
    stop("'legend' must be logical vector of length 1")
  if(legend)
    if(!is.character(desc) || length(desc) != 2)
      stop("'desc' must be a character vector of length 2")
  ## x
  if(!is.numeric(x))
    stop("'x' must be numeric.")
  if(is.matrix(x)){
    if(nrow(x) != length(ind))
      stop("'nrow(x)' must be equal to 'length(ind)'. If you have missing wells, please use the argument 'ind'")
  } else {
    if(length(x) != length(ind))
      stop("'length(x)' must be equal to 'length(ind)'. If you have missing wells, please use the argument 'ind'")
    x = matrix(x, ncol=1)
  }
  
  data <- matrix(NA, ncol=ncol(x), nrow=nrwell)
  ## ind
  if (any(duplicated(ind)) || (max(ind, na.rm = TRUE) > nrwell))
    stop("'ind' must be vector of unique indices for vector 'x'")
  data[ind, ] <- x
  wh <- (1:nrwell)[ind]
  whIn <- seq(along=x)
  ## callArgs
  if(!is.null(callArgs)){
    if(default)
      stop("'callArgs' are not allowed for default plotting function")
    if(!is.data.frame(callArgs) || nrow(callArgs)!=nrow(x))
      stop("'callArgs' must be data frame with same number of rows as 'x'")
  }
    
  ## device & font size ##
  if(!missing(device)) {
    height=((width-0.01*width)/(ncol+1))*(nrow+1)
    args <- list(width = width, height = height)
    if (!missing(file))
      args <- append(args, file)
    if(device %in% c("X11", "x11", "windows", "quartz", "gnome", "GTK"))
      if(names(dev.cur())!="null device")
        dev.off()
    do.call(device, args = args)

    ## units to pixel ##
    xlim = c(0, ncol + 1)
    ylim = c(0, nrow + 1)
    fw = diff(xlim)/0.9
    fh = diff(ylim)/0.9
    u2px = function(x) (x - xlim[1])/fw * width
    u2py = function(y) (y - ylim[1])/fh * height

    ## fontsize
    fontsize = ifelse(device %in% c("png", "jpeg"),
      ceiling(12 * width/900),
      ceiling(12 * width/9))
  
  } else {
    usr <- par("usr") # need this to reinitialize plot
    plot.new() 
    par(usr=usr, plt=usr)
    rg <- rectGrob()  # this is a helper grob to determine vp/device width
    width <- as.numeric(convertHeight(grobWidth(rg), "points"))
    height <- as.numeric(convertHeight(grobHeight(rg), "points"))
    height=min(height, (width-0.01*width)/(ncol+1)*(nrow+1))
    vpFrame <- viewport(width=unit(width-2, "points"), height=unit(height-2, "points"))
    pushViewport(vpFrame)  # this vp makes sure we are plotting the correct aspect ratio
    fontsize = ceiling(12 * width/900)
  }
  
  ## create false colors ##
  if(default){
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
    circcol <- matrix(thepalette[z2icol(data)], ncol=ncol(data), nrow=nrow(data))
    ##deal with NA values
    switch(na.action, zero = {
      circcol[is.na(circcol)] <- thepalette[z2icol(0)]
    }, omit = {
      wh <- which(!is.na(circcol))
      whIn <- which(!is.na(x))
    }, xout = {
      nawell <- which(is.na(circcol))
      sel <- nawell[which((nawell %in% ind))]
      circcol[is.na(circcol)] <- "lightgray"
    }, stop(paste("Invalid value of 'na.action':", na.action)))
  }#end if


  
  ## create grid graphic ##
  vp1 <- viewport(width=0.9, x=0, just="left") #main well vp
  ##plot legend
  if(default && legend){
    vp2 <- viewport(width=0.1, x=0.9, just="left") #main legend vp
    pushViewport(vp2)
    vp3 <- viewport(height=0.85, width=0.8) #legend desc vp
    pushViewport(vp3)
    if (any(!is.na(desc)))
      grid.text(y=c(0.95, 0.05), x=0.1, just="left", desc,
                gp=gpar(fontsize=fontsize, cex=1.4, fontface="bold",
                col=thepalette[c(length(thepalette), 1)]))
    vp4 <- viewport(height=0.8, width=0.1, yscale=c(xrange), xscale=c(0,1),
                    x=0.1, just="left") #legend bar vp
    pushViewport(vp4)
    cr <- colorRampPalette(col)
    nb <- 50
    cols <- thepalette[floor(seq(1, nrcolors, length=nb))]
    for(i in 1:nb)
      grid.rect(y=unit(0+i/50, "npc"), height=unit(1/50, "npc"),
                gp=gpar(fill=cols[i], col=cols[i]), just="top")
    grid.rect()
    at <- signif(seq(xrange[1], xrange[2], length=6)[2:5],2)
    grid.yaxis(at=at, gp=gpar(fontsize=fontsize, cex=1), main=FALSE, label=FALSE)
    grid.text(x=unit(3.5, "native"), y=unit(at, "native"), at, rot=90,
              gp=gpar(fontsize=fontsize, cex=1))
    popViewport(3)
  }#end if
  ##plot wells
  pushViewport(vp1)
  vp5 <- viewport(height=0.9, y=0, just="bottom", xscale=c(0, ncol+1),
                  yscale=c(0, nrow+1)) #outer well vp
  pushViewport(vp5)
  x0 <- ((0:(nrwell - 1))%%ncol)
  y0 <- (nrow - (0:(nrwell - 1))%/%ncol)-1
  xpos <- x0[wh]
  ypos <- y0[wh]
  xdat <- as.matrix(data[wh,])
  vp6 <- viewport(width=unit(1-(1/(ncol+1)), "npc"),
                  height=unit(1-(1/(nrow+1)), "npc"),
                  x=unit(1/(ncol+1), "npc"), y=unit(1/(nrow+1), "npc"),
                  just=c("left","bottom"), xscale=c(0, ncol),
                  yscale=c(0, nrow)) #inner well vp
  pushViewport(vp6)
  grid.rect()
  if(default){ #set args for default plotting function 
    xcol <- as.matrix(circcol[wh,])
    callArgs <- data.frame(col=I(xcol))
  }
  for(i in 1:length(xpos)){
    vptemp <- viewport(height=unit(1, "native"), width=unit(1, "native"),
                       x=unit(xpos[i], "native"), y=unit(ypos[i], "native"),
                       just=c("left", "bottom")) #individual well vp
    pushViewport(vptemp)
    thisArgs <- c(list(data=xdat[i,]), as.list(callArgs[i, ]))
    do.call(gridCall, thisArgs) #call plotting function
    grid.text(x=0.5, y=0.5, info[i], gp=gpar(fontsize=fontsize, cex=1.8))
    if(na.action=="xout" & all(is.na(xdat[i,]))){
      grid.lines(unit(c(0.1, 0.9), "native"), unit(c(0.1,
                 0.9), "native"), gp=gpar(lwd=2, col="darkgray")) 
      grid.lines(unit(c(0.9, 0.1), "native"), unit(c(0.1,
                 0.9), "native"), gp=gpar(lwd=2, col="darkgray"))
    }
    popViewport(1)
  }  
  popViewport(1)
  ##plot well description
  vp7 <- viewport(width=unit(1-(1/(ncol+1)), "npc"), height=unit(1/(nrow+1),
            "npc"), x=unit(1/(ncol+1), "npc"), y=unit(1/(nrow+1), "npc"),
            just=c("left","top"), xscale=c(0, ncol), yscale=c(0, 1)) #well horiz. text vp
  pushViewport(vp7)
  grid.text(x=unit(unique(x0)+0.5, "native"), y=unit(0.9, "native"),
            1:ncol, just="top", gp=gpar(fontsize=fontsize, cex=1.6,
                                  fontface="bold"))
  popViewport(1)
  vp8 <- viewport(width=unit(1/(ncol+1), "npc"), height=unit(1-(1/(nrow+1)),
            "npc"), x=unit(1/(ncol+1), "npc"), y=unit(1/(nrow+1), "npc"),
            just=c("right","bottom"), xscale=c(0, 1), yscale=c(0, nrow)) #well vert. text vp
  pushViewport(vp8)
  grid.text(x=unit(0.9, "native"), y=unit(unique(y0)+0.5, "native"),
            LETTERS[1:nrow], gp=gpar(fontsize=fontsize, cex=1.6, fontface="bold"),
            just="right")
  popViewport(3)
  if (!missing(main)){
    vp9 <- viewport(height=0.1, y=0.9, just="bottom") #well header vp
    pushViewport(vp9)
    grid.text(main, gp=gpar(fontsize=fontsize, cex=1.8, fontface="bold"))
  }
  popViewport(2)
  
  ## return value ##
  res <- list(which=wh)

  ## imageMap coordinates ##   
  if(!missing(device)) {
    dev.off()
    if (device %in% c("png", "jpeg")){
      x0 = 1.5 + (wh - 1)%%ncol
      y0 = 0.1 * diff(ylim) + 0.6 + (wh - 1)%/%ncol
      dx = dy = 0.4
      x1 = u2px(x0 - dx)
      x2 = u2px(x0 + dx)
      y1 = u2py(y0 - dy)
      y2 = u2py(y0 + dy)
      coord = floor(cbind(x1, y1, x2, y2) + 0.5)
      res = append(res, coord=coord, height = args$height, width = args$width)
    }
  }

  return(invisible(res))
}





  
.drawPie <- function(data, col1, col2, col3){
  xpos <- ypos <- 0.5
  r=0.4
  rad <- c(0, cumsum(data)/sum(data)*2)
  nredges <- 180
  colors <- c(col1, col2, col3)
  for(i in 2:length(rad)){
    phi <- seq(rad[i-1] * pi , rad[i] * pi, len=ceiling(nredges*rad[i]))
    x <- c(xpos, r * cos(phi)+xpos, xpos)
    y <- c(ypos, r * sin(phi)+ypos, ypos)
    grid.polygon(x,y, gp=gpar(fill=colors[i-1]))  
  }
}




.drawCircle <- function(data){
  grid.circle(0.5, 0.5, r=max(0.1, data/2-0.1), gp=gpar(fill="red")) 
}


plotPlate <- function (x, nrow = 8, ncol = 12, ind = 1: (ncol*nrow), main,
                       xrange, col, device, file, width, na.action = "zero",
                       desc = as.character(c(NA, NA)), char){
    if (!is.numeric(ncol) || length(ncol) != 1)
        stop("'ncol' must be a numeric vector of length 1")
    if (!is.numeric(nrow) || length(nrow) != 1)
        stop("'nrow' must be a numeric vector of length 1")
    nrwell <- ncol * nrow
    if ((length(ind) != length(x)) | (max(ind, na.rm = TRUE) > nrwell))
      stop("'ind' must be vector of unique indices for vector 'x'")
    y <- rep(NA, nrwell)
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

    xlim = c(-0.5, ncol + 1.5)
    ylim = c(-0.5, nrow + 1.5)
    colbarwid = 0.4
    fw = diff(xlim) + colbarwid
    fh = diff(ylim)
    height <- width/fw * fh
    args <- list(width = width, height = height)
    if (!missing(device)) {
        if (!is.character(device))
            stop("'device' must be a character")
        if (!missing(file))
            args <- append(args, file)
        do.call(device, args = args)
    }
    layout(matrix(1:2, ncol = 2), widths = c(diff(xlim), colbarwid),
        heights = 1)
    u2px = function(x) (x - xlim[1])/fw * width
    u2py = function(y) (y - ylim[1])/fh * height
    par(mai = c(0, 0, 0, 0))
    cex = 1.5
    plot(x = 0, y = 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",
        xaxs = "i", yaxs = "i", xlim = xlim, ylim = ylim)
    graphics::text((1:ncol), 0, paste(1:ncol), adj = c(0.5, 0),
        cex = cex)
    graphics::text(0, (nrow:1), LETTERS[1:nrow], adj = c(0, 0.5),
        cex = cex)
    if (!missing(main))
        graphics::text((ncol + 1)/2, nrow + 1, main, adj = c(0.5,
            1), cex = cex)
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
    radius = 0.45
    xc = radius * cos(seq(0, 2 * pi, len = 73))
    yc = radius * sin(seq(0, 2 * pi, len = 73))
    x0 = 1 + (0:(nrwell - 1))%%ncol
    y0 = nrow - (0:(nrwell - 1))%/%ncol
    switch(na.action, zero = {
        circcol[is.na(circcol)] <- thepalette[z2icol(0)]
        wh <- (1:nrwell)[ind]
    }, omit = {
        wh <- which(!is.na(circcol))
    }, xout = {
        #nawell <- which(which((is.na(circcol)) %in% ind)
        nawell <- which(is.na(circcol))
        sel <- nawell[which((nawell %in% ind))]
        circcol[is.na(circcol)] <- "lightgray"
        wh <- (1:nrwell)[ind]
    }, stop(paste("Invalid value of 'na.action':", na.action)))
    for (i in wh) {
        polygon(x = x0[i] + xc, y = y0[i] + yc, col = circcol[i])
    }
  
    if(na.action=="xout")
        for (i in sel){
	  lines(c(x0[i]-0.45, x0[i]+0.45), c(y0[i]-0.45, y0[i]+0.45),lwd=2,
                col="darkgray") 
	  lines(c(x0[i]-0.45, x0[i]+0.45), c(y0[i]+0.45, y0[i]-0.45), lwd=2,
                col="darkgray")
	}
    if(!missing(char))
       for (i in which(!is.na(char)))
	       	text(x0[i], y0[i], char[i], cex=1.5)
    
    xmin = 0.5
    xmax = ncol + 0.5
    ymin = 0.5
    ymax = nrow + 0.5
    polygon(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax,
        ymax, ymin))
    par(mai = 0.5 * c(1, 0, 1, 0.12))
    image(x = 0, y = icol2z(1:nrcolors), z = matrix(1:nrcolors,
        nrow = 1), col = thepalette, xaxt = "n")
    if (any(!is.na(desc))) {
        mtext(desc[1], side = 3, line = 0.5, font = 2, adj = 0.5,
            col = thepalette[nrcolors])
        mtext(desc[2], side = 1, line = 0.5, font = 2, adj = 0.5,
            col = thepalette[1])
    }
    if (!missing(device) && device %in% c("pdf", "png", "jpeg"))
        dev.off()
    x0 = 1 + (wh - 1)%%ncol
    y0 = 1 + (wh - 1)%/%ncol
    dx = dy = 0.4
    x1 = u2px(x0 - dx)
    x2 = u2px(x0 + dx)
    y1 = u2py(y0 - dy)
    y2 = u2py(y0 + dy)
    if (!missing(ind))
        wh <- which(!is.na(x[ind]))
    res <- list(which = wh, coord = floor(cbind(x1, y1, x2, y2) +
        0.5), height = args$height, width = args$width)
    invisible(res)
}

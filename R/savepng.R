## width and height in pixels
savepng <- function(fn, dir, width=480, asp=1) {
  fn <- paste(fn, ".png", sep="")
  if(!missing(dir))
    fn <- file.path(dir, fn)
  dev.copy(png, filename=fn , width=width, height=width*asp)
  dev.off()
  return(fn)
}

## width and height in inches
savepdf <- function(fn, dir, width=6, asp=1) {
  fn = paste(fn, ".pdf", sep="")
  if(!missing(dir))
    fn <- file.path(dir, fn)
  dev.copy(pdf, file=fn, width=width, height=width*asp)
  dev.off()
  return(fn)
}
saveeps <- function(fn, dir, width=6, asp=1) {
  fn = paste(fn, ".eps", sep="")
  if(!missing(dir))
    fn <- file.path(dir, fn)
  dev.copy(postscript, file=fn, width=width, height=width*asp,
           horizontal = FALSE, onefile = FALSE, paper = "special")
  dev.off()
  return(fn)
}
savetiff <- function(fn, dir, density=150, keepeps=TRUE, ...) {
  epsfn  <- saveeps(fn, dir, ...)
  fn     <- paste(fn, ".tiff", sep="")
  if(!missing(dir))
    fn <- file.path(dir, fn)
  cmd   <- paste("convert", epsfn, "-density", density, "-compress LZW", fn)
  cat(cmd, "\n")
  system(cmd)
  if(!keepeps) file.remove(epsfn)
  return(fn)
}
  

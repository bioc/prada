## width and height in pixels
savepng <- function(fn, dir=".", width=480, asp=1) {
  fn = paste(fn, ".png", sep="")
  dev.copy(png, filename = file.path(dir, fn), width=width, height=width*asp)
  dev.off()
  return(fn)
}
## width and height in inches
savepdf <- function(fn, dir=".", width=6, asp=1) {
  fn = paste(fn, ".pdf", sep="")
  dev.copy(pdf, file = file.path(dir, fn), width=width, height=width*asp)
  dev.off()
  return(fn)
}
saveeps <- function(fn, dir=".", width=6, asp=1) {
  fn = paste(fn, ".eps", sep="")
  dev.copy(postscript, file=file.path(dir, fn), width=width, height=width*asp,
           horizontal = FALSE, onefile = FALSE, paper = "special")
  dev.off()
  return(fn)
}
savetiff <- function(fn, dir = ".", density=150, ...) {
  epsfn <- saveeps(fn, dir=dir, ...)
  fn    <- paste(fn, ".tiff", sep="")
  cmd   <- paste("convert -density", density, epsfn, fn)
  cat(cmd, "\n")
  system(cmd)
  return(fn)
}
  

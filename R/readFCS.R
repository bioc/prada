## For specifications of FACS 3.0 see
## http://www.isac-net.org and the file
## fcs3.html in the doc directory

readFCS <- function(filename)
{
  stopifnot(is.character(filename), length(filename)==1, filename!="")
  con <- file(filename, open="rb")
  
  offsets <- readFCSheader(con)
  txt     <- readFCStext(con, offsets)
  dat     <- readFCSdata(con, offsets, txt)

  close(con)
  
  attr(dat, "text") <- txt
  return(dat)
}

readFCSgetPar <- function(x, pnam) {
  stopifnot(is.character(x), is.character(pnam)) 
  i <- match(pnam, names(x))
  if(any(is.na(i)))
    stop(paste("Parameter(s)", pnam, "not contained in 'x'"))
  return(x[i])
}

readFCSheader <- function(con) {
  seek(con, 0)
  version <- readChar(con, 6)
  stopifnot(version %in% c("FCS2.0", "FCS3.0"))
  
  tmp <- readChar(con, 4)
  stopifnot(tmp=="    ")
  
  offsets   <- numeric(6)
  names(offsets) <- c("textstart", "textend", "datastart", "dataend", "anastart", "anaend")
  for(i in 1:length(offsets))
    offsets[i] <- as.integer(readChar(con, 8))
  
  stopifnot(!any(is.na(offsets)))
  stopifnot(!any(offsets[1:4]==0))
  
  return(offsets)
}

readFCStext <- function(con, offsets, delimiter="\\") {
  seek(con, offsets["textstart"])
  txt <- readChar(con, offsets["textend"]-offsets["textstart"]+1)
  sp  <- strsplit(txt, split=delimiter, fixed=TRUE)[[1]]
  if(sp[1]=="") sp<-sp[-1]
  stopifnot(length(sp)%%2==0)
  rv <- sp[seq(2, length(sp), by=2)]
  names(rv) <- sp[seq(1, length(sp)-1, by=2)]
  return(rv)
}

readFCSdata <- function(con, offsets, txt, endian="big") {
  if (readFCSgetPar(txt, "$BYTEORD") != "4,3,2,1")
    stop(paste("Don't know how to deal with $BYTEORD", readFCSgetPar(txt, "$BYTEORD")))
  
  if (readFCSgetPar(txt, "$DATATYPE") != "I")
    stop(paste("Don't know how to deal with $DATATYPE", readFCSgetPar(txt, "$DATATYPE")))
  
  if (readFCSgetPar(txt, "$MODE") != "L")
    stop(paste("Don't know how to deal with $MODE", readFCSgetPar(txt, "$MODE")))

  nrpar    <- as.integer(readFCSgetPar(txt, "$PAR"))
  range    <- as.integer(readFCSgetPar(txt, paste("$P", 1:nrpar, "R", sep="")))
  bitwidth <- as.integer(readFCSgetPar(txt, paste("$P", 1:nrpar, "B", sep="")))

  if (!all(bitwidth==16))
    stop(paste("Don't know how to deal with the bit widths"))

  seek(con, offsets["datastart"])
  dat <- readBin(con, "integer", n=offsets["dataend"]-offsets["datastart"]+1,
                 size=2, signed=FALSE, endian=endian)
  stopifnot(length(dat)%%nrpar==0)
  dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
  colnames(dat) <- readFCSgetPar(txt, paste("$P", 1:nrpar, "N", sep=""))
  return(dat) 
}


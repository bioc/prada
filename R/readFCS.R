## For specifications of FACS 3.0 see
## http://www.isac-net.org and the file
## fcs3.html in the doc directory

readFCS <- function(filename)
{
  stopifnot(is.character(filename), length(filename)==1, filename!="")
  con <- file(filename, open="rb")
  
  offsets <- readFCSheader(con)
  txt     <- readFCStext(con, offsets)
  mat     <- readFCSdata(con, offsets, txt)

  close(con)
  
  if(as.integer(readFCSgetPar(txt, "$TOT"))!=nrow(mat))
    stop(paste("file", filename, "seems to corrupted."))
  
  return(new("cytoFrame", exprs=mat, description=txt))
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
  if(!version %in% c("FCS2.0", "FCS3.0"))
    stop("This does not seem to be a valid FCS2.0 or FCS3.0 file")
  
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

readFCStext <- function(con, offsets) {
  seek(con, offsets["textstart"])
  txt <- readChar(con, offsets["textend"]-offsets["textstart"]+1)
  delimiter <- substr(txt, 1, 1)
  sp  <- strsplit(substr(txt, 2, nchar(txt)), split=delimiter, fixed=TRUE)[[1]]
  if(length(sp)%%2!=0)
    stop("In readFCStext: unexpected format of the text segment")
  rv <- sp[seq(2, length(sp), by=2)]
  names(rv) <- sp[seq(1, length(sp)-1, by=2)]
  return(rv)
}

readFCSdata <- function(con, offsets, x, endian="big") {
  if (readFCSgetPar(x, "$BYTEORD") != "4,3,2,1")
    stop(paste("Don't know how to deal with $BYTEORD", readFCSgetPar(x, "$BYTEORD")))
  
  if (readFCSgetPar(x, "$DATATYPE") != "I")
    stop(paste("Don't know how to deal with $DATATYPE", readFCSgetPar(x, "$DATATYPE")))
  
  if (readFCSgetPar(x, "$MODE") != "L")
    stop(paste("Don't know how to deal with $MODE", readFCSgetPar(x, "$MODE")))

  nrpar    <- as.integer(readFCSgetPar(x, "$PAR"))
  range    <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "R", sep="")))
  bitwidth <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "B", sep="")))
  bitwidth <- unique(bitwidth)
  if(length(bitwidth)!=1)
    stop("Sorry, I am expecting the bitwidth to be the same for all parameters")

  seek(con, offsets["datastart"])

  size <- bitwidth/8
  if (!size %in% c(2, 4, 8))
    stop(paste("Don't know how to deal with bitwidth", bitwidth))

  dat <- readBin(con, "integer", n = (offsets["dataend"]-offsets["datastart"]+1)/size,
                 size=size, signed=FALSE, endian=endian)

  stopifnot(length(dat)%%nrpar==0)
  dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
  colnames(dat) <- readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))
  return(dat) 
}


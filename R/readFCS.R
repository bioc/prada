## For specifications of FACS 3.0 see
## http://www.isac-net.org and the file
## fcs3.html in the vignettes directory

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

  coffs <- character(6)
  for(i in 1:length(coffs))
    coffs[i] <- readChar(con=con, nchars=8)
  
  ioffs <- as.integer(coffs)
  names(ioffs) <- c("textstart", "textend", "datastart", "dataend", "anastart", "anaend")
  
  stopifnot(all(!is.na(ioffs) | coffs=="        "), !any(ioffs[1:4]==0))
  return(ioffs)
}

readFCStext <- function(con, offsets) {
  seek(con, offsets["textstart"])
  txt <- readChar(con, offsets["textend"]-offsets["textstart"]+1)
  txt <- iconv(txt, "", "latin1", sub="byte")
  delimiter <- substr(txt, 1, 1)
  sp  <- strsplit(substr(txt, 2, nchar(txt)), split=delimiter, fixed=TRUE)[[1]]
  ## if(length(sp)%%2!=0)
  ##  stop("In readFCStext: unexpected format of the text segment")
  rv <- sp[seq(2, length(sp), by=2)]
  names(rv) <- sp[seq(1, length(sp)-1, by=2)]
  return(rv)
}

readFCSdata <- function(con, offsets, x) {
  endian <- switch(readFCSgetPar(x, "$BYTEORD"),
    "4,3,2,1" = "big",
        "2,1" = "big",
        "1,2" = "little",     
    "1,2,3,4" = "little",
    stop(paste("Don't know how to deal with $BYTEORD", readFCSgetPar(x, "$BYTEORD"))))
  

  dattype <- switch(readFCSgetPar(x, "$DATATYPE"),
	"I" = "integer",
	"F" = "numeric",
    	stop(paste("Don't know how to deal with $DATATYPE", readFCSgetPar(x, "$DATATYPE")))) 
  
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
  if (!size %in% c(1, 2, 4, 8))
    stop(paste("Don't know how to deal with bitwidth", bitwidth))

  dat <- readBin(con, dattype, n = (offsets["dataend"]-offsets["datastart"]+1)/size,
                 size=size, signed=FALSE, endian=endian)

  stopifnot(length(dat)%%nrpar==0)
  dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
  colnames(dat) <- readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))
  return(dat) 
}



## This is a wrapper for all flow cytometry import function
## in prada and rflowcyt (which is now defunct). Be specifying the objectModel argument
## the user can choose between the different data representations

read.fcs <- function(filename=NULL, objectModel="prada", ...){
  if(objectModel=="prada"){
    require(prada)
    if(!is.null(filename) && length(filename)==1){
      return(readFCS(filename))
    }else{
      return(readCytoSet(filename, ...))
    }
  }
  ## Give defunct message
  if(objectModel=="FCS"){
      .Defunct(msg="The rflowcyt package is no longer part of Bioconductor. Please consider using package 'flowCore'")
  }
  stop("'objectModel' must be either 'prada' or 'FCS'")
}
  

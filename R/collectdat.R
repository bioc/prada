collectdat = function(thedir, verbose=TRUE) {

  gdat = function(i, j=NULL, aschar=FALSE) {
    res = get("dat", envir=v[[i]])
    if(missing(j)) return(res)
  
    res = res[[j]]
    if(aschar)
      res=as.character(res)
    return(res)
  }

  files  = dir(thedir, pattern = "^dorit.*rda$")

  uselesspattern = "-No="
  nruseless = length(grep(uselesspattern, files))
  if(nruseless>0)
    stop(paste("There are", nruseless, "files that contain", uselesspattern,
               ". Please consider deleting them.\n"))
  
  ## create and load environments
  v = lapply(1:length(files), function(i) new.env())
  if(verbose) cat("Loading", length(files), "files: ")
  for (i in seq(along=files)) {
    if (verbose && (i%%50==0)) cat(i, "")
    load(file.path(thedir, files[i]), envir=v[[i]])
  }
  if(verbose) cat("\n")

  theclass   = sapply(gdat(1), class)
  thecname   = colnames(gdat(1))
  if(verbose) cat("Checking consistency of column names and types.\n")
  for (i in 2:length(v)) {
    if(!all(sapply(gdat(i), class) == theclass))
      stop("Not conformable classes")
    if(!all(colnames(gdat(i)) == thecname))
      stop("Mismatch colnames")
  }

  ## merge
  if(verbose) cat("Merging", length(theclass), "columns:\n")
  dat        = vector(mode="list", length=length(theclass))
  names(dat) = thecname
  for (j in 1:length(theclass)) {
    if(verbose) cat(j, ":", theclass[j], " ", sep="")
    dat[[j]] = unlist(lapply(1:length(v), gdat, j=j, aschar=(theclass[j]=="factor")))
  }
  if(verbose) cat("\n")

  ## expRepeat: character -> integer
  if(verbose) cat("Converting columns to integers.\n")
  dat[["expRepeat"]] = as.all(dat[["expRepeat"]], "integer")

  ## numeric -> integer
  for (i in c("wellSub", "Area (pixels)", "Perimeter", "Projection x", "Projection y",
              "Longest Segment Length", "x1Left", "y1Top", "x2Right", "y2Bottom",
              "Number of Holes", "Hole's Area (pixels)", "Field")) {
    stopifnot(i %in% names(dat))
    dat[[i]] <- as.all(dat[[i]], "integer")
  }

  ## repair the plateSub field
  cat("plateSub as read from XML files:")
  j1 = which(names(dat)=="plateSub") ## new
  j2 = which(names(dat)=="wellSub")  ## new
  stopifnot(length(j1)==1, length(j2)==1)

  print(table(dat[[j1]]))
  s = gsub("c-1-2|c_1-2", "c12", as.character(dat[[j1]]))
  s = gsub("c-3-4|c_3-4", "c34", s)
  s = gsub("_1-48",  "1-48", s)
  s = gsub("_49-96", "49-96", s)
  cat("plateSub repaired:")
  print(table(s))

  ## merge plateSub, wellSub --> well
  well = mapSub2Plate(s, dat[[j2]]) 
  dat[[j1]] = well
  names(dat)[j1] = "well"
  dat = dat[-j2]   ## remove j2

  ## rename column headers
  nm <- names(dat)
  newnam <- c("brdu",      "trsf",              "dapi")
  oldnam <- c("Brdu mean", "transfection mean", "Dapi mean")
  for (i in seq(along=newnam)) {
    j <- which(nm==oldnam[i])
    stopifnot(length(j)==1)
    names(dat)[j] <- newnam[i]
  }
  
  ## convert to data frame and save
  cat("Converting to data frame.\n")
  dat = data.frame(dat)

  file <- file.path(thedir, paste("dat", gsub(" ", "-", date()), "rda", sep="."))
  cat("Saving", file, "\n")
  save(dat, file=file)
  
  ## check
  l1        <- sapply(1:length(v), function(i) { nrow(get("dat", envir=v[[i]])) } )
  names(l1) <- sapply(1:length(v), function(i) {
    x <- get("dat", envir=v[[i]])
    x <- paste(as.character(x$expId), as.character(x$expRepeat), as.character(x$plateSub),
               as.character(x$wellSub), as.character(x$cloneId))
    x <- unique(x)
    stopifnot(length(x)==1)
    return(x)
  } )
  l2 = table(paste(as.character(dat$expId), as.character(dat$expRepeat),
                   as.character(dat$well), as.character(dat$cloneId)))
  stopifnot(sort(l1)==sort(l2))
}  
  

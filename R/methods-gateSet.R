## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
  signature="gateSet", definition=function(object) {
    validObject(object)
    if(object@name == "ALL" && all(names(object) == "ALL")){
      cat("Initial gateSet object containing all observations\n")
    }else{
      gnames <- sapply(object@glist, names, simplify=TRUE)
      gtypes <- glogics <- gcnames <- NULL
      if(length(object@glist)>0){
        for(i in 1:length(object@glist)){
          gtypes <- c(gtypes, object@glist[[i]]@type)
          glogics <- c(glogics, object@glist[[i]]@logic)
          gcnames <- c(gcnames, paste(object@glist[[i]]@colnames,
                                      collapse=" and "))
        }
        msg <- paste("gateSet object '", object@name, "' containing ",
                     length(object@glist), " individual gate(s):\n", sep="")
        tab <- cbind("Name:"=gnames, "Type:"=gtypes, "Logic:"=glogics,
                     "Variables:"=gcnames)
        rownames(tab) <- as.character(1:length(gnames))
        cat(msg)
        print(format(tab, just="left"))
      }
    }
    return(NULL)}, valueClass="character")
## ==========================================================================

## ==========================================================================
## names method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
  signature="gateSet", definition=function(x) {
    names <- NULL
    for(i in 1:length(x))
      names <- c(names, x@glist[[i]]@name)
    return(names)}, valueClass="character")
## ==========================================================================

## ==========================================================================  
## repacement method for names
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("names",
  signature=c("gateSet"), definition=function(x, value) {
    if(!is.character(value) || length(value) != length(x))
      stop("\nreplacement attribute must be same length as object")
    for(i in 1:length(x))
      x@glist[[i]]@name <- value[i]
    validObject(x)
    return(x)})
## ==========================================================================

## ==========================================================================
## length method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("length",
  signature="gateSet", definition=function(x) {
    return(length(x@glist))}, valueClass="integer")
## ==========================================================================

## ==========================================================================
## subsetting method to gate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("[[",
  signature="gateSet",
  definition=function(x, i, j="missing", drop="missing") {
    if(!missing(j))
      stop("invalid number of dimensions")
    if(length(i)!=1)
      stop("Subsetting to single items only")
    if(!i %in% 1:length(x@glist))
      stop("Subset out of bounds")
    ret <- x@glist[[i]]
    validObject(ret)
    return(ret)
   },
   valueClass="gate")
## ==========================================================================

## ==========================================================================
## subsetting method to gateSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("[",
  signature="gateSet",
  definition=function(x, i, j="missing", drop="missing") {
    if(!missing(j))
      stop("invalid number of dimensions")
    if(!all(i %in% 1:length(x@glist)))
      stop("Subset out of bounds")
    x@glist <- x@glist[i]
    validObject(x)
    return(x)
   },
   valueClass="gateSet")
## ==========================================================================

## ==========================================================================
## append gates method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("appendGates",
          signature=c("gateSet"),
          definition=function(x, ...){
            arglist <- list(...)
            if(!all(sapply(arglist, is, "gate")))
              stop("Can only append 'gate' objects")
            gsnames <- names(x@glist)
            gnames <- sapply(arglist, names)
            x@glist <- c(x@glist, arglist)
            names(x@glist) <- c(gsnames, gnames)
            validObject(x)
            return(x)},
          valueClass="gateSet")
  

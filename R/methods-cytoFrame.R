## ==========================================================================
## accessor method for slot exprs
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("exprs",
  signature="cytoFrame", definition=function(object) object@exprs,
  valueClass="matrix")
## ==========================================================================


## ==========================================================================
## replace method for slot exprs
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("exprs",
  signature=c("cytoFrame", "matrix"), definition=function(object, value) {
    object@exprs <- value
    return(object)})
## ==========================================================================


## ==========================================================================
## accessor method for slot description
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("description",
  signature="cytoFrame", definition=function(object) object@description,
  valueClass="character")
## ==========================================================================


## ==========================================================================
## replace method for slot description
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("description",
  signature=c("cytoFrame", "character"), definition=function(object, value) {
    object@description <- value
    return(object)})
## ==========================================================================


## ==========================================================================
## accessor method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames",
  signature="cytoFrame", definition=function(x, do.NULL="missing",
  prefix="missing") colnames(exprs(x)), valueClass="character")
## ==========================================================================


## ==========================================================================
## replace method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("colnames",
  signature=c("cytoFrame", "ANY"), definition=function(x, value) {
    colnames(x@exprs) <- value
    return(x)})
## ==========================================================================


## ==========================================================================
## plot method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("plot",
  signature(x="cytoFrame", y="missing"),
  definition=function(x, col=densCols(exprs(x)[,1:2]), pch=20, ...){
    values=exprs(x)
    plot(values, col=col, pch=pch, ...)})
## ==========================================================================


## ==========================================================================
## the $-operator
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"$.cytoFrame" <- function(x, val)
    (description(x))[val]
## ==========================================================================

## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
  signature="cytoFrame", definition=function(object) {
    dm <- dim(exprs(object))
    msg <- paste("cytoFrame object with ", dm[1], " cells and ", dm[2],
    " observables:\n", paste(colnames(exprs(object)), collapse=" "), 
    "\nslot 'description' has ", length(description(object)),
    " elements\n", sep="")
    if(length(object@gate)!=1 || (object@gate@name!="ALL" &&
       object@gate@glist[[1]]@type!="ALL"))
      msg <- paste(msg, "gate(s) assigned to frame in gateSet '",
                   object@gate@name, "':\n   ",
                   paste("'", names(object@gate),
                   collapse="'  ", sep=""), "'\n", sep="")
    cat(msg)            
    return(msg)}, valueClass="character")
## ==========================================================================


## ==========================================================================
## subsetting method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("[",
  signature="cytoFrame", definition=function(x, i, j, ..., drop=FALSE) {
    exprs(x) <-  switch(1+missing(i)+2*missing(j),
         { exprs(x)[i, j, ..., drop=drop] },
         { exprs(x)[ , j, ..., drop=drop] },
         { exprs(x)[i,  , ..., drop=drop] },
         { exprs(x)[ ,  , ..., drop=drop] } )
    x
  },
  valueClass="cytoFrame")
## ==========================================================================


## ==========================================================================
## accessor method for slot gate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("gate",
  signature=signature("cytoFrame"),
    definition=function(x) {
    return(x@gate)})
## ==========================================================================


## ========================================================================== 
## replace method for slot gate with gate objects
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("gate", signature=c("cytoFrame", "gate"),
  definition=function(object, value) {
    glist <- list(value)
    names(glist) <- names(value)
    sname=names(value)
    gset <- new("gateSet", name=sname, glist=glist)
    gate(object) <- gset
    validObject(object)
    return(object)})      
## ==========================================================================


## ========================================================================== 
## replace method for slot gate with gateSet objects
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("gate", signature=c("cytoFrame", "gateSet"),
  definition=function(object, value){
    gnames <- names(value)
    names(value@glist) <- gnames
    object@gate <-  value
    return(object)})      
## ==========================================================================


## ==========================================================================
## drawGate method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("drawGate",
  signature=signature("cytoFrame"),
    definition=function(x, ...) {
    g <- gateMatrix(exprs(x), ...)
    return(invisible(g))})
## ==========================================================================


## ==========================================================================
## nrow method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("nrow",
  signature=signature("cytoFrame"),
    definition=function(x) {
    return(nrow(x@exprs))})
## ==========================================================================


## ==========================================================================
## ncol method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("ncol",
  signature=signature("cytoFrame"),
    definition=function(x) {
    return(ncol(x@exprs))})
## ==========================================================================


## ==========================================================================
# applyGate method on matrix with gate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("matrix", "gate"),
  definition=function(data, x) {
    return(data[x@gateFun(data),])}, valueClass="matrix")
## ==========================================================================


## ==========================================================================
## applyGate method on matrix with gateSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("matrix", "gateSet"),
  definition=function(data, x) {
    ## test for validity of objects
    if(!is.matrix(data) || (!is.double(data) && !is.integer(data)))
      stop("'data' must be real or integer matrix")
    cnData <- colnames(data)  
    gfuns <- glogics  <- NULL
    cnGates <- list() 
    for(i in 1:length(x)){
      gfuns <- c(gfuns, x@glist[[i]]@gateFun)
      glogics <- c(glogics, x@glist[[i]]@logic)
      cnGates[[i]] <- x@glist[[i]]@colnames  
    }
    cng <- unique(unlist(cnGates))
    miss <- which(!cng %in% cnData) 
    if(length(miss)!=0)
      stop("Need variable(s) '", paste(cng, collapse="' and '"),
           "' in colnames of data matrix to apply gate\n'",
           paste(cnData[miss], collapse="' and '"),
           "' not present", sep="")
    if(length(cng)<1)
       return(invisible(!logical(nrow(data))))
 
    fCalls <- paste("gfuns[[", 1:length(gfuns), "]](data[,",
                    paste("c(\"", lapply(cnGates, paste,
                    collapse="\", \""), "\")", sep=""), ", drop=FALSE])",
                    sep="")
    or <- which(glogics=="|")
    if(length(or)!=0)
      andCalls <- fCalls[-c(or, or-1)]
    else
      andCalls <- fCalls
    orCalls <- fCalls[sort(unique(c(or, or-1)))]
    andSel <- eval(parse(text=paste(andCalls, collapse=" & ")))
    orSel <-  eval(parse(text=paste(orCalls, collapse=" | ")))
    if(length(andSel) & length(orSel)){
      allSel <- andSel & orSel
    }else if(length(andSel)){
      allSel <- andSel
    }else{
      allSel <- orSel
    }
    return(data[allSel,, drop=FALSE])}, valueClass="matrix")
## ==========================================================================


# applyGate method on cytoFrame with gate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("cytoFrame", "gate"),
  definition=function(data, x) {
    exprs(data) <- applyGate(exprs(data), x)
    return(data)
  }, valueClass="cytoFrame")
## ==========================================================================


# applyGate method on cytoFrame with gateSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("cytoFrame", "gateSet"),
  definition=function(data, x) {
    exprs(data) <- applyGate(exprs(data), x)
    return(data)
  })
## ==========================================================================


## ==========================================================================
## applyGate method on cytoFrame with character
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("cytoFrame", "character"),
  definition=function(data, x) {
    if(!all(x %in% names(data@gate)))
      stop("\ngate not assigned to this object\n  available gates: '",
           paste(names(data@gate), collapse="' '", sep=""), "'")
    wh <- match(x, names(data@gate))
    gate <- new("gateSet", glist=data@gate@glist[wh], name="tmp")
    
    exprs(data) <- applyGate(exprs(data), gate)
    return(data)
  }, valueClass="cytoFrame")
## ==========================================================================


## ==========================================================================
## applyGate method on cytoFrame with logical
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("cytoFrame", "logical"),
  definition=function(data, x) {
    if(!x)
      return(data)
    gate <- data@gate
    exprs(data) <- applyGate(exprs(data), gate)
    return(data)
  }, valueClass="cytoFrame")
## ==========================================================================


## ==========================================================================
## applyGate method on cytoFrame with numeric
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("cytoFrame", "numeric"),
  definition=function(data, x) {
    x <- as.integer(x)
    if(!all(x %in% 1:length(data@gate)))
      stop("gate index out of bounds")
    gate <- data@gate[x]
    exprs(data) <- applyGate(exprs(data), gate)
    return(data)
  }, valueClass="cytoFrame")
## ==========================================================================


## ==========================================================================
## range method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
range.cytoFrame <- function(..., na.rm=FALSE){
  arglist=list(...)
  if(length(arglist)!=1)
    stop("too many arguments")
  range(exprs(arglist[[1]]), na.rm=na.rm)
}
## ==========================================================================

## ==========================================================================
## append gates method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("appendGates",
          signature=c("cytoFrame"),
          definition=function(x, ...){
            arglist <- list(...)
            if(!all(sapply(arglist, is, "gate")))
              stop("Can only append 'gate' objects")
            gsnames <- names(x@gate@glist)
            gnames <- sapply(arglist, names)
            x@gate@glist <- c(x@gate@glist, arglist)
            names(x@gate@glist) <- c(gsnames, gnames)
            validObject(x)
            return(x)},
          valueClass="cytoFrame")
## ==========================================================================

## ==========================================================================
## hist method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("hist",
  signature(x="cytoFrame"),
  definition=function(x, ...){
    x=exprs(x)
    hist(x, ...)})
## ==========================================================================

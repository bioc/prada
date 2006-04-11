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
    gset <- new("gateSet", name=names(value), glist=glist)
    gate(object) <- gset
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
## applyGate method on cytoFrame with character
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("character", "cytoFrame"),
  definition=function(x, data) {
    if(!all(x %in% names(data@gate)))
      stop("\ngate not assigned to this object\n  available gates: '",
           paste(names(data@gate), collapse="' '", sep=""), "'")
    wh <- match(x, names(data@gate))
    gate <- new("gateSet", glist=data@gate@glist[wh], name="tmp")
    
    exprs(data) <- applyGate(gate, exprs(data))
    return(data)
  }, valueClass="cytoFrame")
## ==========================================================================

## ==========================================================================
## applyGate method on cytoFrame with logical
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("logical", "cytoFrame"),
  definition=function(x, data) {
    if(!x)
      return(data)
    gate <- data@gate
    exprs(data) <- applyGate(gate, exprs(data))
    return(data)
  }, valueClass="cytoFrame")
## ==========================================================================

## ==========================================================================
## applyGate method on cytoFrame with numeric
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("numeric", "cytoFrame"),
  definition=function(x, data) {
    x <- as.integer(x)
    if(!all(x %in% 1:length(data@gate)))
      stop("gate index out of bounds")
    gate <- data@gate[x]
    exprs(data) <- applyGate(gate, exprs(data))
    return(data)
  }, valueClass="cytoFrame")
## ==========================================================================


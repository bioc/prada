## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
  signature="gate", definition=function(object) {
    validObject(object)
    if(object@name=="ALL" && object@type=="ALL"){
      msg <-  paste("Initial gate object including all observations\n")
    }else{
      logic <- switch(object@logic, "&"="AND",
                      "|"="OR",
                      stop("unknown gate logic!"))
      msg <- paste("gate object '", object@name, "' of type '", object@type,
                   "' applied to variables '",
                   paste(object@colnames, collapse="' and '"),
                   "'\n  logic combination with other gates: ", logic, "\n",
                   sep="", collapse="")
    }
    cat(msg)
    return(msg)}, valueClass="character")
## ==========================================================================


## ==========================================================================  
## names method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
  signature="gate", definition=function(x) {
    return(x@name)}, valueClass="character")
## ==========================================================================


## ==========================================================================  
## repacement method for names
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("names",
  signature=c("gate"), definition=function(x, value) {
    if(!is.character(value) || length(value) != 1)
      stop("\nreplacement attribute must be same length as object")
      x@name <- value
    validObject(x)
    return(x)})
## ==========================================================================


# as.gateSet method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("as.gateSet",
  signature=signature("gate"),
  definition=function(from) {
    theClass <- class(from)
    if(!(theClass=="gate"))
      stop("Can't coerce object of class '", theClass, "' to class 'gateSet'")
    gl <- list(from)
    names(gl) <- from@name
    ret <- new("gateSet", glist=gl, name=from@name)
    validObject(ret)
    return(ret)
  }, valueClass="gateSet")
## ==========================================================================


## ==========================================================================  
## lines method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("lines",
  signature="gate", definition=function(x, ...) {
    if(!x@type %in% c("polygon", "rectangle"))
      stop("Don't know how to deal with gate of type '", x@type, "'")
    if(nrow(x@boundaries)>0)
      lines(x@boundaries, ...)})
## ==========================================================================

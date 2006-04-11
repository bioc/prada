## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
  signature="gate", definition=function(object) {
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


# applyGate method on matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("gate", "matrix"),
  definition=function(x, data) {
    return(data[x@gateFun(data),])}, valueClass="matrix")
## ==========================================================================  


# applyGate method on cytoFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("gate", "cytoFrame"),
  definition=function(x, data) {
    exprs(data) <- applyGate(x, exprs(data))
    return(data)
  }, valueClass="cytoFrame")
## ==========================================================================


# as.gate method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("as.gateSet",
  signature=signature("gate"),
  definition=function(x) {
    ret <- new("gateSet", glist=list(x), name=x@name)
    return(ret)
  }, valueClass="gateSet")
## ==========================================================================

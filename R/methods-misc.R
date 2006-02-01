## ==========================================================================
## drawGate method for matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("drawGate",
  signature=signature("matrix"),
          definition=function(x, ...) {
    g <- gateMatrix(x, ...)
    return(g)})
## ==========================================================================

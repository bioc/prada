## ==========================================================================
## drawGate method for matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("drawGate",
  signature=signature("matrix", "ANY"),
          definition=function(x, ...) {
    g <- gateMatrix(x, ...)
    return(g)})
## ==========================================================================

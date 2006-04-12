## ==========================================================================
## append gates to gateSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
combineGates <- function(gate, ..., name="gateSet1"){
  arglist <- c(list(gate), list(...))
  if(!all(sapply(arglist, is, "gate")))
     stop("Can only combine 'gate' objects")
  names(arglist) <- sapply(arglist, names)
  return(new("gateSet", name=name, glist=arglist))
}

  

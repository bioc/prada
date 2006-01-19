## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
  signature="gateSet", definition=function(object) {
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
## applyGate method on matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("gateSet", "matrix"),
  definition=function(x, data) {
    ## test for validity of objects
    if(!is.matrix(data) || (!is.real(data) && !is.integer(data)))
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
    if(length(cng)<3)
       return(invisible(!logical(nrow(data))))
 
    fCalls <- paste("gfuns[[", 1:length(gfuns), "]](data[,",
                    paste("c(\"", lapply(cnGates, paste,
                    collapse="\", \""), "\")", sep=""), "])", sep="")
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
    return(invisible(allSel))}, valueClass="logical")
## ==========================================================================

# applyGate method on cytoFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyGate",
  signature=signature("gateSet", "cytoFrame"),
  definition=function(x, data) {
    return(applyGate(x, exprs(data)))
  })
## ==========================================================================

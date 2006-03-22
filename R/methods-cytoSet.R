## ==========================================================================
## accessor method for slot phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("phenoData",
  signature="cytoSet", definition=function(object) object@phenoData,
  valueClass="phenoData")
## ==========================================================================

## ========================================================================== 
## replace method for slot phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("phenoData", c("cytoSet", "phenoData"),
  function(object, value) {
    object@phenoData <- value
    return(object)})      
## ==========================================================================

## ==========================================================================
##  accessor method for slot phenoData@pData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("pData",
  signature="cytoSet", definition=function(object) object@phenoData@pData,
  valueClass="data.frame")
## ==========================================================================

## ==========================================================================
## accessor method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames",
  signature="cytoSet",
  definition=function(x, do.NULL="missing", prefix="missing") x@colnames,
  valueClass="character")
## ==========================================================================

## ==========================================================================
## replace method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("colnames",
  signature=c("cytoSet", "ANY"), definition=function(x, value) {
    x@colnames <- value
    return(x)})
## ==========================================================================

## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
  signature="cytoSet",
  definition=function(object) {
  cat("cytoSet object with", length(object), "cytoFrames and",
      "colnames\n", paste(colnames(object)), "\n")})
## ==========================================================================

## ==========================================================================
## length method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("length",signature(x="cytoSet"),
          function(x) nrow(x@phenoData@pData))
## ==========================================================================

## ==========================================================================
## subsetting method to cytoSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("[",
  signature="cytoSet",
  definition=function(x, i, j="missing", drop="missing") {
    fr <-new.env(hash=TRUE)
    if(is.numeric(i)||is.logical(i)) {
      nm <- x@phenoData$name[i]
    } else {
      nm <- i
      i <- match(nm, x@phenoData$name)
    }
    multiassign(nm, mget(nm, x@frames), envir=fr, inherits=FALSE) 
    new("cytoSet",
      frames=fr,
      phenoData=x@phenoData[i, ],
      colnames=colnames(x))
   },
   valueClass="cytoSet")
## ==========================================================================

## ==========================================================================
## subsetting method to cytoFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("[[",
  signature="cytoSet",
  definition=function(x, i, j="missing") {
    if(length(i)!=1)
      stop("subscript out of bounds (index must have length 1 in '[[')")
    if(is.numeric(i))
      i <- x@phenoData$name[i]
    rv <- get(i, x@frames, inherits=FALSE)
    colnames(exprs(rv)) <- x@colnames
    return(rv)
  },
  valueClass="cytoFrame")
## ==========================================================================

## ==========================================================================
## replace method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("[[", signature="cytoSet",
  definition=function(x, i, j, ..., value) {
    if(!is.matrix(value))
      stop("'value' must be a matrix.")
    if(ncol(value)!=length(x@colnames))
      stop(paste("'value' has", ncol(value),
                 "columns, while the other elements of this cytoSet have",
                 length(x@colnames), "- these numbers should be the same."))
    if(!is.null(colnames(value)))
      if(!all(colnames(value)==x@colnames))
        stop("'value' must have the same colnames as the other elements",
             "of this cytoSet.")
    if(length(i)!=1)
      stop("subscript out of bounds (index must have length 1 in '[[<-')")
    if(is.numeric(i))
      i <- x@phenoData$name[i]
    theFrame <- get(i, envir=x@frames)
    exprs(theFrame) <- value
    colnames(theFrame) <- NULL
    assign(i, theFrame, envir=x@frames, inherits=FALSE)
    return(x)
  },
  valueClass="cytoFrame")
## ==========================================================================

## ==========================================================================
## plot method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("plot",
  signature(x="cytoSet", y="missing"),
  definition=function(x, col, main, gate, pch=20, ...){
    dcol <- missing(col)
    dm <- missing(main)
    sel <- TRUE
    frames <- format(pData(x)[,"name", drop=FALSE])
    cat("Available frames in cytSet:\n")
    print(frames)
    userAnswer <- NULL
    while(is.null(userAnswer)){
      userAnswer <- readline("Plot which frame? ('a' for all): ")
      if(userAnswer == "a"){
          for(i in 1:length(x)){
            values <- exprs(x[[i]])
            if(!missing(gate))
              sel <- applyGate(gate, values)
            if(dcol)
              col <- densCols(values[sel,1:2])
            if(dm)
              main <- paste("frame #", i, " (", frames[i,], ")", sep="") 
            plot(values[sel,], main=main, col=col, pch=pch, ...)
            par(ask=TRUE)
          } #end for
       }else{
         if(! userAnswer %in% as.character(1:length(x))){
           userAnswer <- NULL
           cat("Invalid entry!")
         }else{
           values=exprs(x[[as.integer(userAnswer)]])
           if(!missing(gate))
              sel <- applyGate(gate, values)
           if(dcol)
             col=densCols(values[sel,1:2])
           if(dm)
             main <- paste("frame #", userAnswer, " (",
                           frames[as.integer(userAnswer),], ")", sep="") 
           plot(values[sel,], col=col, main=main, pch=pch, ...)
         } #end else
       } #end else
    } #end while
    par(ask=FALSE)
  })
## ==========================================================================






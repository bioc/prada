
#################### Class definitions ####################################
#cytoFrame
setClass("cytoFrame",                
  representation(exprs="matrix",
                 description="character",
                 well="integer"),
  prototype=list(exprs=matrix(numeric(0), nrow=0, ncol=0),
                 description=c(note="empty"), well=as.integer(1)),
 validity=function(object){
   is.matrix(object@exprs)&&is.character(object@description)&&
   is.integer(object@well)&&length(object@well)==1
 })


#cytoSet
setClass("cytoSet",                   
  representation(frames="environment",
                 phenoData="phenoData",
                 colnames="character"),
  prototype=list(frames=new.env(),
                 phenoData=new("phenoData",
                   pData=data.frame(name=I(character(0))),
                   varLabels=list(name="Name in frame")),
                 colnames=character(0)),
  validity=function(object){
    nc <- length(colnames(object))
    is(object@phenoData, "phenoData") &&
    is(object@colnames, "character") &&
    is(object@frames, "environment") &&
    "name" %in% colnames(pData(object@phenoData)) &&
    setequal(ls(object@frames, all.names=TRUE), object@phenoData$name) &&
    all(sapply(ls(object@frames, all.names=TRUE), function(x)
      { fr <- get(x, envir=object@frames, inherits=FALSE)
       is(fr, "cytoFrame") && is.null(colnames(fr))  &&
        ncol(exprs(fr))==nc } ))
  })
############################################################################




############################### Generic functions ##########################
## replace colnames
if(!isGeneric("colnames<-"))
  setGeneric("colnames<-", function(x, value)
    standardGeneric("colnames<-"))

## set colnames
if(!isGeneric("colnames"))
  setGeneric("colnames", function(x, do.NULL=TRUE, prefix="col")
    standardGeneric("colnames"))

## plot
if(!isGeneric("plot"))
  setGeneric("plot", useAsDefault=plot)

##These are already defined as generic functions in Biobase:
##if(!isGeneric("exprs<-"))
##  setGeneric("exprs<-", function(object, value)
##    standardGeneric("exprs<-"))
##
##if(!isGeneric("exprs"))
##  setGeneric("exprs", function(object)
##    standardGeneric("exprs"))
##
##if(!isGeneric("description<-"))
##  setGeneric("description<-", function(object, value)
##    standardGeneric("description<-"))
##
##if(!isGeneric("description"))
##  setGeneric("description", function(object)
##    standardGeneric("description"))
##
##if(!isGeneric("length"))
##  setGeneric("length", function(x)
##   standardGeneric("colnames"))

######################################################################



####################### formal methods ###############################
## accessor method for slot exprs
setMethod("exprs",
  signature="cytoFrame",
  definition=function(object) object@exprs,
  valueClass="matrix")

## accessor method for slot description
setMethod("description",
  signature="cytoFrame",
  definition=function(object) object@description,
  valueClass="character")

## accessor method for slot colnames cytoFrame
setMethod("colnames",
  signature="cytoFrame",
  definition=function(x, do.NULL="missing", prefix="missing")
          colnames(exprs(x)),
  valueClass="character")

## accessor method for slot colnames cytoSet
setMethod("colnames",
  signature="cytoSet",
  definition=function(x, do.NULL="missing", prefix="missing") x@colnames,
  valueClass="character")

## accessor method for slot phenoData
setMethod("phenoData",
  signature="cytoSet",
  definition=function(object) object@phenoData,
  valueClass="phenoData")

##  accessor method for slot phenoData@pData
setMethod("pData",
  signature="cytoSet",
  definition=function(object) object@phenoData@pData,
  valueClass="data.frame")

## replace method for slot exprs
setReplaceMethod("exprs",
   signature=c("cytoFrame", "matrix"),
   definition=function(object, value) {
     object@exprs <- value
     return(object)
   })

## replace method for slot description
setReplaceMethod("description",
   signature=c("cytoFrame", "character"),
   definition=function(object, value) {
     object@description <- value
     return(object)
   })

## replace method for slot colnames cytoFrame
setReplaceMethod("colnames",
  signature=c("cytoFrame", "ANY"),
  definition=function(x, value) {
    colnames(x@exprs) <- value
    return(x)
  })

## replace method for slot colnames cytoSet
setReplaceMethod("colnames",
  signature=c("cytoSet", "ANY"),
  definition=function(x, value) {
    x@colnames <- value
    return(x)
  })

## replace method for slot phenoData
setReplaceMethod("phenoData", c("cytoSet", "phenoData"),
  function(object, value) {
    object@phenoData <- value
    return(object)
  })      

## plot method for cytoFrames
setMethod("plot",
  signature(x="cytoFrame", y="missing"),
  definition=function(x, col=densCols(exprs(x)[,1:2]), pch=20, ...){
    values=exprs(x)
    plot(values, col=col, pch=pch, ...)
  })

## the $-operator
"$.cytoFrame" <- function(x, val)
    (description(x))[val]

## show method for cytoFrame
setMethod("show",
  signature="cytoFrame",
  definition=function(object) {
    dm <- dim(exprs(object))
    msg <- paste("cytoFrame object with ", dm[1], " cells and ", dm[2],
       " observables:\n", paste(colnames(exprs(object)), collapse=" "), 
       "\nslot 'description' has ", length(description(object)),
       " elements\n", sep="")
    cat(msg)
    return(msg)
  },
  valueClass="character")

## show method for cytoSet
setMethod("show",
  signature="cytoSet",
  definition=function(object) {
  cat("cytoSet object with", length(object), "cytoFrames and",
      "colnames\n", paste(colnames(object)), "\n")
})

## length method for cytoSet
setMethod("length",signature(x="cytoSet"),
          function(x) nrow(x@phenoData@pData))

## subsetting method for cytoFrame
setMethod("[",
  signature="cytoFrame",
  definition=function(x, i, j, ..., drop=FALSE) {
    exprs(x) <-  switch(1+missing(i)+2*missing(j),
         { exprs(x)[i, j, ..., drop=drop] },
         { exprs(x)[ , j, ..., drop=drop] },
         { exprs(x)[i,  , ..., drop=drop] },
         { exprs(x)[ ,  , ..., drop=drop] } )
    x
  },
  valueClass="cytoFrame")


## subsetting method for cytoSet 2 cytoSet
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

## subsetting method for cytoSet 2 cytoFrame
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

## replace method for cytoSet
setReplaceMethod("[[", signature="cytoSet",
  definition=function(x, i, j, ..., value) {
    if(!is.matrix(value))
      stop("'value' must be a matrix.")
    if(ncol(value)!=length(x@colnames))
      stop(paste("'value' has", ncol(value), "columns, while the other elements of this cytoSet have",
                 length(x@colnames), "- these numbers should be the same."))
    if(!is.null(colnames(value)))
      if(!all(colnames(value)==x@colnames))
        stop("'value' must have the same colnames as the other elements of this cytoSet.")
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
################################################################################







library(Biobase)

setClass("cytoFrame",
  representation(exprs="matrix",
                 description="character"),
  prototype=list(exprs=matrix(numeric(0), nrow=0, ncol=0),
                 description=c(note="empty")),
 validity=function(object){
   is.matrix(object@exprs)&&is.character(object@description)
 })

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

## accessor methods for slots exprs and description
setMethod("exprs",      "cytoFrame", function(object) object@exprs)
setMethod("description", "cytoFrame", function(object) object@description)

setReplaceMethod("exprs", c("cytoFrame", "matrix"),
   function(object, value) {
     object@exprs <- value
     return(object)
   })

setReplaceMethod("description", c("cytoFrame", "character"),
   function(object, value) {
     object@description <- value
     return(object)
   })

"$.cytoFrame" <- function(x, val)
    (description(x))[val]

setMethod("show", "cytoFrame", function(object) {
    dm <- dim(exprs(object))
    cat("cytoFrame object with", dm[1], "cells and", dm[2], "observables:\n")
    cat(colnames(exprs(object)), "\n")
    cat("slot 'description' has length", length(description(object)), "\n")
})

setMethod("[", "cytoFrame", function(x, i, j, ..., drop=FALSE) {
  exprs(x) <-  switch(1+missing(i)+2*missing(j),
         { exprs(x)[i, j, ..., drop=drop] },
         { exprs(x)[ , j, ..., drop=drop] },
         { exprs(x)[i,  , ..., drop=drop] },
         { exprs(x)[ ,  , ..., drop=drop] } )
  x
})

setReplaceMethod("[", "cytoFrame", function(x, i, j, ..., value) {
  exprs(x)[i, j, ...] <- value
  x
  })


setClass("cytoSet",
  representation(frames="environment",
                 phenoData="phenoData",
                 colnames="character"),
  prototype=list(frames=new.env(),
                 phenoData=new("phenoData",
                   pData=data.frame(framename=I(character(0))),
                   varLabels=list(framename="Name in frame")),
                 colnames=character(0)),
  validity=function(object){
    is(object@phenoData, "phenoData") &&
    is(object@colnames, "character") &&
    is(object@frame, "environment") &&
    "framename" %in% colnames(pData(object@phenoData)) &&
    setequal(ls(object@frame), object@phenoData$wellname) &&
    all(sapply(ls(object@frame), function(x)
      { is(get(x, inherits=FALSE), "cytoFrame") &&
        colnames(exprs(get(x, inherits=FALSE)))==NULL })) &&
    all(sapply(ls(object@frame), ncol)==length(object@colnames))
  })

setMethod("[", "cytoSet", function(x, i, j, ..., drop=FALSE) {
  if(!missing(j))
    stop("incorrect number of dimensions")
  fr <-new.env(hash=TRUE)
  nm <- object@phenoData$framename[i]
  multiassign(nm, mget(nm, object@frames), envir=fr, inherits=FALSE) 
  new("cytoSet",
      frames=fr,
      phenoData=x@phenoData[i, ],
      colnames=x@phenoData)
}, valueClass="cytoSet")

setMethod("[[", "cytoSet", function(x, i, j, ...) {
  if(!missing(j))
    stop("incorrect number of dimensions")
  if(length(i)!=1)
    stop("subscript out of bounds (index must have length 1 in '[[')")
  rv <- get(x@phenoData$framename[i], x@frames, inherits=FALSE)
  colnames(exprs(rv)) <- x@colnames
  return(rv)
}, valueClass="cytoFrame")

## show method for cytoSet
setMethod("show", "cytoSet", function(object) {
  cat("cytoSet object which contains", nrow(object@phenoData), "cytoFrames ",
      "with colnames\n", colnames(object), "\n")
  show(phenoData(object))
})

## get and set colnames slot of cytoSet
if(!isGeneric("setColnames<-"))
  setGeneric("setColnames<-", function(object, value)
    standardGeneric("setColnames<-"))

if(!isGeneric("getColnames"))
  setGeneric("getColnames", function(object)
    standardGeneric("getColnames"))

setMethod("getColnames", "cytoSet", function(object) object@colnames)
setReplaceMethod("setColnames", c("cytoSet", "character"),
  function(object, value) {
    object@colnames <- value
    return(object)
  })      

## get and set phenoData slot of cytoSet
setMethod("phenoData", "cytoSet", function(object) object@phenoData)
setReplaceMethod("phenoData", c("cytoSet", "phenoData"),
  function(object, value) {
    object@phenoData <- value
    return(object)
  })      

          

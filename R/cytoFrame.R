library(Biobase)

setClass("cytoFrame",
  representation(exprs="matrix",
                 description="character"),
  prototype=list(exprs=matrix(numeric(0), nrow=0, ncol=0),
                 description=c(note="empty")),
 validity=function(object){
   is.matrix(object@exprs)&&is.character(object@description)
 })

## Generic functions
if(!isGeneric("colnames<-"))
  setGeneric("colnames<-", function(x, value)
    standardGeneric("colnames<-"))

if(!isGeneric("colnames"))
  setGeneric("colnames", function(x, do.NULL=TRUE, prefix="col")
    standardGeneric("colnames"))

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
setMethod("exprs",
  signature="cytoFrame",
  definition=function(object) object@exprs,
  valueClass="matrix")

setMethod("description",
  signature="cytoFrame",
  definition=function(object) object@description,
  valueClass="character")

setMethod("colnames",
  signature="cytoFrame",
  definition=function(x, do.NULL="missing", prefix="missing")
          colnames(exprs(x)),
  valueClass="character")

## replace methods
setReplaceMethod("exprs",
   signature=c("cytoFrame", "matrix"),
   definition=function(object, value) {
     object@exprs <- value
     return(object)
   })

setReplaceMethod("description",
   signature=c("cytoFrame", "character"),
   definition=function(object, value) {
     object@description <- value
     return(object)
   })

setReplaceMethod("colnames",
  signature=c("cytoFrame", "ANY"),
  definition=function(x, value) {
    colnames(x@exprs) <- value
    return(x)
  })

## the $-operator
"$.cytoFrame" <- function(x, val)
    (description(x))[val]

## show method
setMethod("show",
  signature="cytoFrame",
  definition=function(object) {
    dm <- dim(exprs(object))
    msg <- paste("cytoFrame object with", dm[1], "cells and", dm[2],
       "observables:\n", colnames(exprs(object)), "\n",
       "slot 'description' has length", length(description(object)),
       "\n", sep="")
    cat(msg)
    return(msg)
  },
  valueClass="character")

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

## FIXME: do we need this or is it odd?
##
## setReplaceMethod("[", "cytoFrame", function(x, i, j, ..., value) {
##   exprs(x)[i, j, ...] <- value
##  x
##  })

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
    nc <- length(colnames(object))
    TRUE||
    is(object@phenoData, "phenoData") &&
    is(object@colnames, "character") &&
    is(object@frames, "environment") &&
    "framename" %in% colnames(pData(object@phenoData)) &&
    setequal(ls(object@frames), object@phenoData$framename) &&
    all(sapply(ls(object@frames), function(x)
      { fr <- get(x, envir=object@frames, inherits=FALSE)
        is(fr, "cytoFrame") && is.null(colnames(fr))  &&
        ncol(exprs(fr))==nc } ))
  })

setMethod("[",
  signature="cytoSet",
  definition=function(x, i, j="missing", drop="missing") {
    fr <-new.env(hash=TRUE)
    nm <- x@phenoData$framename[i]
    multiassign(nm, mget(nm, x@frames), envir=fr, inherits=FALSE) 
    new("cytoSet",
      frames=fr,
      phenoData=x@phenoData[i, ],
      colnames=x@phenoData)
   },
   valueClass="cytoSet")

setMethod("[[",
  signature="cytoSet",
  function(x, i, j="missing") {
    if(length(i)!=1)
      stop("subscript out of bounds (index must have length 1 in '[[')")
    rv <- get(x@phenoData$framename[i], x@frames, inherits=FALSE)
    colnames(exprs(rv)) <- x@colnames
    return(rv)
  },
  valueClass="cytoFrame")

## show method for cytoSet
setMethod("show",
  signature="cytoSet",
  definition=function(object) {
  cat("\tcytoSet object. Its colnames are:\n\t",
      colnames(object), "\n", sep="")
  show(phenoData(object))
})

setMethod("colnames",
  signature="cytoSet",
  definition=function(x, do.NULL="missing", prefix="missing") x@colnames,
  valueClass="character")

setReplaceMethod("colnames",
  signature=c("cytoSet", "ANY"),
  definition=function(x, value) {
    x@colnames <- value
    return(x)
  })      

## get and set phenoData slot of cytoSet
setMethod("phenoData", "cytoSet", function(object) object@phenoData)
setReplaceMethod("phenoData", c("cytoSet", "phenoData"),
  function(object, value) {
    object@phenoData <- value
    return(object)
  })      

          

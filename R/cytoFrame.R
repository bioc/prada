setClass("cytoFrame",
         representation(exprs="matrix",
                        description="character"))

## accessor methods for slots exprs and description
if(!isGeneric("exprs<-"))
  setGeneric("exprs<-", function(object, value)
    standardGeneric("exprs<-"))

if(!isGeneric("exprs"))
  setGeneric("exprs", function(object)
    standardGeneric("exprs"))

if(!isGeneric("description<-"))
  setGeneric("description<-", function(object, value)
    standardGeneric("description<-"))

if(!isGeneric("description"))
  setGeneric("description", function(object)
    standardGeneric("description"))
   
setMethod("exprs",      "cytoFrame", function(object) object@exprs)
setMethod("description", "cytoFrame", function(object) object@description)

setReplaceMethod("exprs", c("cytoFrame", "matrix"),
   function(object, value) {
     object@exprs <- value
     object })

setReplaceMethod("description", c("cytoFrame", "character"),
   function(object, value) {
     object@description <- value
     object })

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

## Biobase
## if (!isClass("exprList"))
##    setClassUnion("exprList", c("list", "environment"))


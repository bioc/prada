.onLoad <- function(lib, pkg){
  require(methods)
  require(Biobase)

  ##These are already defined as generic functions in Biobase:
  if( !isGeneric("description") && !exists("description", mode="function"))  
    setGeneric("description",     function(object) standardGeneric("description"))
  
  if( !isGeneric("description<-") && !exists("description<-", mode="function"))
    setGeneric("description<-",   function(object, value) standardGeneric("description<-"))

  if( !isGeneric("exprs") && !exists("exprs", mode="function"))  
    setGeneric("exprs",           function(object) standardGeneric("exprs"))
  
  if( !isGeneric("exprs<-") && !exists("exprs<-", mode="function"))
    setGeneric("exprs<-",         function(object, value) standardGeneric("exprs<-"))
  
  if( !isGeneric("phenoData") && !exists("phenoData", mode="function"))  
   setGeneric("phenoData",       function(object) standardGeneric("phenoData"))
  
  if( !isGeneric("phenoData<-") && !exists("phenoData<-", mode="function")) 
    setGeneric("phenoData<-",     function(object, value) standardGeneric("phenoData<-"))

  if( !isGeneric("pData") && !exists("pData", mode="function"))    
    setGeneric("pData",           function(object) standardGeneric("pData"))


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
  
## replace method for slot phenoData
setReplaceMethod("phenoData", c("cytoSet", "phenoData"),
  function(object, value) {
    object@phenoData <- value
    return(object)
  })      

  

}

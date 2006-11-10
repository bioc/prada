require("Biobase")


#### FIXME: Might be useful to create a virtual gate class with subclasses
#### for different types of gates (polygon, elliptical, rectangular...)

## ===========================================================================
## gate
## ---------------------------------------------------------------------------
## An object describing a selection applied to a data matrix. Consist of
## a functions that return logical vectors subsetting the data
## ---------------------------------------------------------------------------
setClass("gate",
  representation(name="character",
                 gateFun="function",
                 colnames="character",
                 logic="character",
                 type="character",
                 boundaries="matrix" ## this is important for plotting
                 ),
  prototype=list(name="ALL", gateFun=function(x) TRUE,
                 logic="&", type="unknown", boundaries=matrix(ncol=2, nrow=0)),
  validity=function(object){
    msg <- TRUE
    if(!is.character(object@colnames) ||
       length(object@colnames)<1)
      msg <- "\nslot 'colnames' must be character vector longer then 1"
    test <- matrix(1:length(object@colnames), ncol=length(object@colnames))
    colnames(test) <- object@colnames
    if(!is.function(object@gateFun) ||
       !is.logical(object@gateFun(test)))
      msg <- paste("\nslot 'gateFun' must be function returning logical",
                   "vector when applied to data matrix")
    if(!is.character(object@name) ||
       length(object@name)!=1)
      msg <- "\nslot 'name' must be character vector of length 1"
    if(!is.character(object@logic) ||
       length(object@logic)!=1 ||
       !object@logic %in% c("&", "|"))
      msg <- "\nslot 'logic' must be character, either '&' or '|'"
    if(!is.character(object@type) ||
       length(object@type)!=1)
      msg <- "\nslot 'type' must be character of length 1"
    return(msg)
  })

## ===========================================================================
## gateSet
## ---------------------------------------------------------------------------
## An object describing a set of individual gating functions to subset
## data, possibly in several dimensions.
## ---------------------------------------------------------------------------
setClass("gateSet",      
  representation(name="character",
                 glist="list"),
  prototype=list(name="ALL", glist=list(ALL=new("gate"))),
  validity=function(object){
    msg <- TRUE
    if(! is.character(object@name) ||
       length(object@name)!=1)
        msg <- "\nslot 'name' must be character vector of length 1"
    if(!is.list(object@glist) ||
        !all(sapply(object@glist, is, "gate")))
        msg <- paste("\nslot 'glist' must be list of length > 0",
                    "with items of class 'gate'")
    gnames <- sapply(object@glist, names)
    if(length(unique(gnames))!=length(gnames))
      msg <- "names of individual gates are not unique"
    if(is.null(names(object@glist)) ||
       length(setdiff(names(object@glist), gnames)) != 0)
      msg <- "names of 'glist' must match names of individual gates"
    return(msg)
  })
        

## ===========================================================================
## cytoFrame
## ---------------------------------------------------------------------------
## A container for flow cytometry measurements with slots exprs, description
## and well. Exprs contains measurement values, description contains 
## information from file headers of FCS file and well contains well position
## on microtiter plate from experiment. Gate contains an object of class
## gateSet, which may be assessed for subsequent operations, e.g.plotting.
## ---------------------------------------------------------------------------
setClass("cytoFrame",                
  representation(exprs="matrix",
                 description="character",
                 well="integer",
                 gate="gateSet"),
  prototype=list(exprs=matrix(numeric(0), nrow=0, ncol=0),
                 description=c(note="empty"), well=as.integer(1),
                 gate=new("gateSet")),
 validity=function(object){
   msg <- TRUE
   if(!is.matrix(object@exprs))
      msg <- "\nslot 'exprs' must be matrix"
   if(!is.character(object@description))
      msg <- "\nslot 'description' must be character vector"
   if(! is.integer(object@well) ||
      length(object@well)!=1)
     msg <- "\nslot 'description' must be integer vector of length 1"
   if(!is(object@gate, "gateSet"))
     msg <- "\nslot 'gate' must be object of class gateSet"
 })

## ===========================================================================
## cytoSet
## ---------------------------------------------------------------------------
## A collection of several cytoFrames making up one experiment. Slots 
## frames, phenoData, colnames. Frames contains the cytoFrame objects,
## phenoData the experiment meta data and colnames the channel names. 
## ---------------------------------------------------------------------------
setClass("cytoSet",                   
  representation(frames="environment",
                 phenoData="AnnotatedDataFrame",
                 colnames="character"),
  prototype=list(frames=new.env(),
                 phenoData=new.env(),
                 phenoData=new("AnnotatedDataFrame",
                   data=data.frame(name=I(character(0))),
                   varMetadata=data.frame(labelDescription="Name in frame", row.names="name")),
                 colnames=character(0)),
  validity=function(object){
    nc <- length(colnames(object))
    is(object@phenoData, "AnnotatedDataFrame") &&
    is(object@colnames, "character") &&
    is(object@frames, "environment") &&
    "name" %in% colnames(pData(object@phenoData)) &&
    setequal(ls(object@frames, all.names=TRUE), object@phenoData$name) &&
    all(sapply(ls(object@frames, all.names=TRUE), function(x)
      { fr <- get(x, envir=object@frames, inherits=FALSE)
       is(fr, "cytoFrame") && is.null(colnames(fr))  &&
        ncol(exprs(fr))==nc } ))
  })

## ===========================================================================
## ddCtSet
## ---------------------------------------------------------------------------
## A subclass to the virtual eSet class to store rtPCR data that has been
## analyzed using the ddCt method.
## ---------------------------------------------------------------------------
setClass("ddCtSet", contains = "eSet",
         validity=function(object){
           msg <- TRUE
           msg <- assayDataValidMembers(assayData(object), c("exprs", "level.err",
                                        "Ct", "Ct.error", "dCt", "dCt.error", "ddCt",
                                        "ddCt.error", "Difference", "numberNA",
                                        "number", "Plate"))
           if(!all(apply(sapply(assayData(object), dim), 2,
                         function(x) identical(dim(assayData(object)[[1]]), x))))
             msg <- "all items of assayData must have same dimesions" 
           return(msg)})

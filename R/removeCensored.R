.removeCensored <- function(x, values, columns, na.rm=TRUE) {
  if(missing(columns))
    columns <- 1:ncol(x)
  if(is.character(columns))
    columns <- match(columns, colnames(x))
  if(!is.numeric(columns) || any(is.na(columns)) ||
     any(columns<1) || any(columns>ncol(x)))
    stop("One or more elements of 'columns' not found in 'x'")

  if(missing(values))
    values <- range(x)

  throwOut <- rep(FALSE, nrow(x))
  for(v in values)
    throwOut <- throwOut | (rowSums(x[, columns, drop=FALSE] == v) > 0)
  if(na.rm)
    throwOut <- throwOut | (rowSums(is.na(x[, columns, drop=FALSE])) > 0)
  
  return(x[!throwOut,, drop=FALSE])
}

if(!isGeneric("removeCensored"))
  setGeneric("removeCensored",  function(x, values, columns, na.rm=TRUE)
    standardGeneric("removeCensored"))

setMethod("removeCensored",
  signature="matrix", 
  valueClass="matrix", 
  definition=.removeCensored)

setMethod("removeCensored",
  signature="data.frame", 
  valueClass="data.frame", 
  definition=.removeCensored)

setMethod("removeCensored",
  signature="cytoFrame",
  valueClass="cytoFrame",
  definition=function(x, values, columns, na.rm=TRUE) {
    exprs(x) <- .removeCensored(exprs(x), values, columns, na.rm)
    x
})


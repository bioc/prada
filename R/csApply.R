csApply <- function(X, FUN, ..., simplify = TRUE)
{
  if (!is(X, "cytoSet"))
    stop("'X' must be of class cytoSet")
  sapply(X@phenoData$name, function(i){
         Y <- exprs(X[[i]])
         attr(Y, "well") <- X[[i]]@well
         FUN(Y, ...)}, simplify=simplify)
}


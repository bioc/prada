csApply <- function (X, FUN, ...)
{
  if (!is(X, "cytoSet"))
    stop("'X' must be of class cytoSet")
  sapply(X@phenoData$name, function(i)
           FUN(exprs(X[[i]]), ...))
}


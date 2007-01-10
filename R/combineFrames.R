combineFrames <- function(x, by) {
  if(!is(x, "cytoSet"))
    stop("'x' must be a cytoSet.")
  if(!is.factor(by))
    by<-factor(by)
  if(length(by)!=length(x))
    stop("Length of 'by' must be the same as length of 'x'.")
  e  <- new.env(hash=TRUE)
  df <- data.frame(name=I(levels(by)))
  vl <- data.frame(labelDescription=I(c(name="Grouping factor, obtained from argument 'by' of 'combineFrames'")))
  for(i in 1:ncol(pData(x))) {
    s <- split(pData(x)[[i]], by)
    s <- lapply(s, unique)
    if(all(listLen(s)==1)) {
      df <- cbind(df, unlist(s))
      vl <- rbind(vl, I(varMetadata(phenoData(x))[i,]), deparse.level=0)
      rownames(vl) <- c(rownames(vl)[-(nrow(vl))],
                        varLabels(phenoData(x))[[i]]) 
      colnames(df)[ncol(df)] <- names(vl)[length(vl)] <- colnames(pData(x))[i]
    }
  }
  for(lev in df$name) {
    cf <- new("cytoFrame", exprs=do.call("rbind", args=csApply(x[which(by==lev)], function(z) z, simplify=FALSE)))
    colnames(cf) <- NULL
    assign(lev, cf, envir=e)
  }
  #rownames(vl) <- varLabels(phenoData(x))
  new("cytoSet", frames=e,
      phenoData=new("AnnotatedDataFrame", data=df, varMetadata=vl),
      colnames=colnames(x))
}



getDye <- function(x) {
  stopifnot(is.character(x))
  cfp1 <- grep("(^C)|(^[A-Za-z][A-Za-z]C)", x)  ## C at the beginning
  cfp2 <- grep("C[0-9]*\\+*[A-Za-z]*$", x)      ## 'C' at the end
  
  yfp1 <- grep("(^Y)|(^[A-Za-z][A-Za-z]Y)", x)  ## Y at the beginning
  yfp2 <- grep("Y[0-9]*\\+*[A-Za-z]*$", x)      ## Y at the end

  yfp2 <- setdiff(yfp2, cfp1) ## C at beginning overrides Y at the end
  cfp2 <- setdiff(cfp2, yfp1) ## Y at beginning overrides C at the end
  
  yfp  <- union(yfp1, yfp2)
  cfp  <- union(cfp1, cfp2)

  stopifnot(length(intersect(cfp, yfp))==0)

  res <- character(length(x))
  res[is.na(x)] <- as.character(NA)
  res[cfp] <- "cfp"
  res[yfp] <- "yfp"
  stopifnot(all(res %in% c("cfp", "yfp") | is.na(res)))
  return(res)  
}

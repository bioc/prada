getPradaPar <- function(parname) {
  stopifnot(is.character(parname), length(parname)==1)
  w <- which(parname==names(getFromNamespace("pradaPars", ns="prada")))
  if (length(w)!=1)
     stop(paste("The parameter", parname, "is not defined"))
  return(getFromNamespace("pradaPars", ns="prada")[[w]])
}

setPradaPars <- function(pars) {
  stopifnot(is.list(pars), is.character(names(pars)))
  assignInNamespace("pradaPars", pars, ns="prada")
  invisible(NULL)
}

pradaPars <- as.list(NULL)

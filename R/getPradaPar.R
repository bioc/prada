getPradaPar <- function(parname) {
  stopifnot(is.character(parname), length(parname)==1)
  penv <- as.environment("package:prada")
  if(!exists("pradaPars", envir=penv))
    stop("Global options 'pradaPars' are not defined.")
  w <- which(parname==names(get("pradaPars", envir=penv)))
  if (length(w)!=1)
     stop(paste("The parameter", parname, "is not defined"))
  return(get("pradaPars", envir=penv)[[w]])
}

setPradaPars <- function(pars) {
  stopifnot(is.list(pars), is.character(names(pars)))
  assign("pradaPars", pars, envir=as.environment("package:prada"))
  invisible(NULL)
}

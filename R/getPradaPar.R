getPradaPar <- function(parname) {
  stopifnot(is.character(parname), length(parname)==1)
  if(!exists("pradaPars", envir=globalenv()))
    stop("Global options 'pradaPars' not defined in the global environment")
  w <- which(parname==names(get("pradaPars", envir=globalenv())))
  if (length(w)!=1)
     stop(paste("The parameter", parname, "is not defined"))
  get("pradaPars", envir=globalenv())[[w]]
}

setPradaPars <- function(pars) {
  stopifnot(is.list(pars), is.character(names(pars)))
  assign("pradaPars", pars, envir=globalenv())
  invisible(NULL)
}

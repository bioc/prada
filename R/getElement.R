getElement = function(node, name, type="numeric") {
  stopifnot(class(node)=="XMLNode")
  nn = names(node)
  stopifnot(length(nn)==2 && all(nn == c("Name", "Val")))
  stopifnot(is.character(name))
  stopifnot(is.character(type))
  
  ch = xmlChildren(node)
  stopifnot(xmlValue(ch[[1]]) == name)
  res = xmlValue(ch[[2]])
  if (type=="character") {
    stopifnot(xmlName(node) == "String")
  } else if (type=="numeric") {
    stopifnot(xmlName(node) %in% c("SGL", "DBL", "I32"))
    res = as.numeric(gsub(",", ".", res))
  } else {
    stop(paste("type=", type, " not known.", sep=""))
  }
  return(res)
}

getVector = function(node, name, elementnames, elementtypes) {
  stopifnot(class(node)=="XMLNode")
  stopifnot(is.character(name) || is.null(name))
  
  n = length(elementnames)
  if(length(elementtypes)==1)
    elementtypes = rep(elementtypes, n)

  if(length(elementtypes)!=n)
    stop(paste("length(elementtypes)=", length(elementtypes), " != ",
               "length(elementnames)=", n, sep=""))
  
  stopifnot(xmlName(node) == "Cluster")
  ch = xmlChildren(node)
  if(!is.null(name))
    stopifnot(xmlValue(ch[[1]]) == name)
  stopifnot(xmlValue(ch[[2]]) == n)
  
  z = vector(mode="list", length=n)
  for (i in 1:n) {
    typ = elementtypes[[i]]
    stopifnot(length(typ)>=1)
    if (length(typ)==1) {
      if (! typ %in% c("character", "numeric"))
        stop(paste("elementtype[[", i, "]] should be one of",
             "\"character\", \"numeric\", but is", typ))
      ## atomic object
      z[[i]] = getElement(ch[[2+i]], name=elementnames[[i]], type=elementtypes[[i]])
    } else {
      ## descend into the nesting
      z[[i]] = getVector(ch[[2+i]], name=NULL, elementnames=elementnames[[i]], elementtypes=elementtypes[[i]])
    }
  }

  z = unlist(z)
  names(z) = unlist(elementnames)
  return(z)
}



splitWellName = function(x) {
  stopifnot(is.character(x) && length(x)==1)
  
  s = strsplit(x, "--")[[1]]
  if(length(s)!=2)
    stop(paste("Malformed wellname ", x, ": expecting only one occurence of --", sep=""))
    
  re = regexpr("\\_", s[1])
  rv = c(substr(s[1], start=1,    stop=re-1),
         substr(s[1], start=re+1, stop=re+1),
         substr(s[1], start=re+2, stop=nchar(s[1])),
         sub("-.*", "", s[2]),  ## wellSub
         sub(".*-", "", s[2]))  ## cloneId

  names(rv) = c("expId", "expRepeat", "plateSub", "wellSub", "cloneId")
  return(rv)
}




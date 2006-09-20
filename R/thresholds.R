thresholds = function(x, y, xthr, ythr){
  if(!missing(y)){
    xn <- cbind(x, y)
  }else{
    xn <- x
  }
  if(!is.numeric(xthr) | !is.numeric(ythr))
    stop("'xthr' and 'ythr' must be numeric vectors")
 
  bx  <- (xn[,1]<=xthr)
  by  <- (xn[,2]<=ythr)
  nbx <- !bx
  nby <- !by
  return(rbind(c(sum( bx&nby), sum(nbx&nby)),
               c(sum( bx& by), sum(nbx& by))))
} 


 

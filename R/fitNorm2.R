fitNorm2 <- function(x, y=NA, scalefac=1, method="covMcd") {

  if (!(is.matrix(x) && ncol(x)==2)){
    if (!length(x)==length(y) || !is.numeric(x) || !is.numeric(y))
      stop("'x' and 'y' must be numeric vectors of equal length")
    x <- cbind(x, y)
  }
  if (nrow(x)<50)
    stop("Not enough data points for reliable analysis")
  if (!is.numeric(scalefac))
    stop("'scalefac' must be numeric")

  cov <- switch(method,
    covMcd = {
      nmax <- 50000
      if (nrow(x)>nmax)
        covMcd(x[sample(nrow(x), nmax),])
      else
        covMcd(x)
    },
    cov.rob = {
      cov.rob(x)
    },
    stop("'method' must be one of 'covMcd' or 'cov.rob'")
 ) ## end of switch
  
 mu   <- cov$center
 S    <- cov$cov
 Sinv <- solve(S)
 w    <- rbind(x[,1], x[,2])-mu
 z    <- Sinv %*% w
 p    <- exp(-0.5 * (z[1,]*w[1,] +  z[2,]*w[2,]))
 sel  <- p > exp(-0.5 * scalefac^2)

 return(invisible(list(mu=mu, S=S, p=p, sel=sel, scalefac=scalefac, data=x)))
}


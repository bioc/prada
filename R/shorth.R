shorth <- function(z) {
  sz    <- sort(z)
  width <- round(0.5*length(z))
  diffs <- sz[(width+1):length(z)] - sz[1:(length(z)-width)]
  q     <- which(diffs==min(diffs))

  # deal with ties: if they lie within 5%, simply take the average
  # otherwise, generate an error 
  if (max(q)-min(q) <= 0.05*length(z)) {
    q <- mean(q)
  } else {
    stop(paste("Error in shorth.r: found more than one q: ", q, "\n", sz[q], "\n", sz[q+width], "\n"))
  }
  return(mean(sz[q:(q+width-1)]))
}


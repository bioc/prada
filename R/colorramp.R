colorramp = function (col)  {
  coord <- as.data.frame(t(col2rgb(col))/255)
  x <- seq(0, 1, length = length(col))
  r <- approxfun(x, coord$red)
  g <- approxfun(x, coord$green)
  b <- approxfun(x, coord$blue)
  function(n) {
    x <- seq(0, 1, length = n)
    rgb(r(x), g(x), b(x))
  }
}

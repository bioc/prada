openHTMLpage = function(name, title="") {
  con = file(paste(name, ".html", sep=""), open="wt")
  writeLines(paste("<html><head><title>", title, "</title></head><body>", sep=""), con)
  return(con)
}

closeHTMLpage = function(con) {
  writeLines("</body></html>", con)
  close(con)
}


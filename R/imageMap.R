imageMap = function(con, imgname, coord, tooltips, url, target="extra") {
  stopifnot(is.character(imgname) && length(imgname)==1)
  stopifnot(is.matrix(coord) && ncol(coord)==4)
  n <- nrow(coord)
  stopifnot(is.character(tooltips) && length(tooltips)==n)
  stopifnot(is.character(url) && length(url)==n)
  stopifnot(is.character(target), !any(is.na(target)))
  stopifnot(!any(is.na(coord)), !any(is.na(tooltips)), !any(is.na(url)))
  
  mapname <- paste("map", gsub(" |/|#", "_", imgname), sep="_")
  base::writeLines(paste("<IMG SRC=\"", imgname, "\" USEMAP=\#", mapname, " BORDER=0>",
                         sep=""), con)
  base::writeLines(paste("<MAP NAME=\"", mapname, "\">", sep=""), con)
  base::writeLines(paste("<AREA SHAPE=\"rect\" HREF=\"", url, "\" TITLE=\"", tooltips,
                   "\" COORDS=\"", coord[,1], ",", coord[,2], ",", coord[,3], ",",
                   coord[,4], "\" TARGET=\"", target, "\">", sep=""), con)
  base::writeLines("</MAP>", con)
}


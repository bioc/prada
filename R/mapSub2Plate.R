mapSub2Plate = function(sub, well) {
  stopifnot(length(sub)==length(well))
  stopifnot(is.character(sub) & is.integer(well))
  stopifnot(all(is.na(well)==is.na(sub)))
  
  is.chsl  = (sub %in% c("a", "b", "c12", "c34"))
  is.other = (sub %in% c("1-48", "49-96", ""))
  if(!all(is.chsl | is.other))
    stop("Undefined values in 'sub'")

  ## chamberslide
  ir = (well-1) %%4     ## 0..3
  ic = (well-1) %/% 4   ## 0..7
  ir = ifelse((ic%%2) > 0, 3-ir, ir)
  sel = (sub == "b");    ir[sel] <- ir[sel]+4
  sel = (sub == "c12");  ic[sel] <- ic[sel]+8
  sel = (sub == "c34");  ic[sel] <- ic[sel]+4; ir[sel] <- ir[sel]+4
  well[is.chsl] = (12*ir+ic+1)[is.chsl]

  ## for other, just pass through

  well <- as.integer(well)
  stopifnot(all(well>=1 & well<=96))
  return(well)
}


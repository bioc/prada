plotPlateDB = function(platename, channel, opt) {

  sqlcmd = paste("SELECT id FROM t_plate_analysis WHERE platename=", inquotes(platename), ";")
  res = sqlQuery(channel, sqlcmd)
  assert.data.frame(res)
  plateid = res$id
  
  sqlcmd = paste("SELECT platewellnr, noofcells, tfefficiency, effect, logp, plateid FROM t_well_analysis",
                 "WHERE plateid = ", plateid, "ORDER BY platewellnr;")
  res = sqlQuery(channel, sqlcmd, as.is = c(4))
  assert.data.frame(res, onerow=FALSE)

  signs = (match(res$effect, c("inh", "no", "act")) - 2)
  stopifnot(all(signs %in% c(-1,0,1)))
  res$values = pmin(-res$logp, 10) * signs
               
  asp = plot.plate(list(z      = pmin(-res$logp, 10) * signs,
                        zlim   = c(-10, 10),
                        donuts = NULL,
                        wellnr = res$platewellnr), opt=opt, name = platename)
  
  png.effect = savepng(paste(platename, "effect", sep="_"), opt$outdir, width=768, asp=asp)
  pdf.effect = savepdf(paste(platename, "effect", sep="_"), opt$outdir, asp=asp)

  asp = plot.plate(list(z    = res$noofcells,
                        zlim = c(0, max(res$noofcells)),
                        donuts = 0.8 * sqrt(1-res$tfefficiency),
                        wellnr = res$platewellnr), opt=opt, name = platename)

  png.noofcells = savepng(paste(platename, "noofcells", sep="_"), opt$outdir, width=768, asp=asp)
  pdf.noofcells = savepdf(paste(platename, "noofcells", sep="_"), opt$outdir, asp=asp)

  sqlcmd = paste("UPDATE t_plate_analysis SET fn_effect=", inquotes(png.effect),
                                           ", fn_noofcells=", inquotes(png.noofcells),
                 "WHERE plateid=", plateid, ";")
  res = sqlQuery(channel, sqlcmd)
  
}  

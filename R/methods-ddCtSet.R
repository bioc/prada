setMethod("initialize", "ddCtSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   exprs = new("matrix"),
                   level.err = new("matrix"),
                   Ct = new("matrix"),
                   Ct.error = new("matrix"),
                   dCt = new("matrix"),
                   dCt.error = new("matrix"),
                   ddCt = new("matrix"),
                   ddCt.error = new("matrix"),
                   Difference = new("matrix"),
                   numberNA = new("matrix"),
                   number = new("matrix"),
                   Plate = new("matrix"),
                   ...)
                   {
            callNextMethod(.Object,
                           assayData = assayDataNew(
                             exprs = exprs,
                             level.err = level.err,
                             Ct = Ct,
                             Ct.error = Ct.error,
                             dCt = dCt,
                             dCt.error = dCt.error,
                             ddCt = ddCt,
                             ddCt.error = ddCt.error,
                             Difference = Difference,
                             numberNA = numberNA,
                             number = number,
                             Plate = Plate,
                             ..., storage.mode="list"
                            ),
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
          })


setMethod("exprs", "ddCtSet", function(object) assayDataElement(object,"exprs"))

setReplaceMethod("exprs", c("ddCtSet","matrix"),
                 function(object, value){
                   if(!identical(dim(value), dim(assayDataElement(object,"exprs"))))
                     stop("replacement value has wrong dimensions")
                   assayDataElementReplace(object, "exprs", value)
                   validObject(object)
                   return(object)
                   })

ddCt <- function(raw.table,
                 calibrationSample,
                 housekeepingGene,
                 sampleInformation=NULL,
                 filename="warning.output.txt"){
  withCallingHandlers({
    if (! all(c("Ct","Sample","Detector","Platename") %in% colnames(raw.table)))
      stop ("Your .txt file must include columns with the following names : 'Ct','Sample','Detector','Platename'.")
  
    aaa <- raw.table$Ct
    bbb <- raw.table$Platename
    reduced.set <- raw.table[,c("Sample","Detector")]

    if (! housekeepingGene %in% reduced.set[,2])
      stop(paste("Your reference gene",housekeepingGene,"is not part of your your .txt file."))
    if (! all(calibrationSample %in% reduced.set[,1]))
      stop("At least one of your reference samples is not part of your .txt file.")

    the.difference <- function(x) return (ifelse (any(is.na(x)), NA, max(diff(sort(x)))))
    sum.na <- function(x) {sum(is.na(x))}
    unique.plate <- function(x) {
      if (length(unique(x))!=1) warning(paste("g-s comb. on more than one plate:",paste(unique(x),collapse=",")))
      return (unique(x)[1])}
    
    number.of.na          <- tapply(aaa,reduced.set,sum.na)               # Nr. der nicht ausgewerteten Punkte
    the.Ct.values         <- tapply(aaa,reduced.set,median,na.rm=TRUE)    # der Median der Triplets
    the.mad.values        <- tapply(aaa,reduced.set,mad,na.rm=TRUE,con=1) # der MAD in einem Triplet
    the.difference.values <- tapply(aaa,reduced.set,the.difference)       # Verhältnis lange zu kurze Strecke
    the.plate             <- tapply(bbb,reduced.set,unique.plate)
    
    for (sample in calibrationSample)
      for( Detector in unique(reduced.set[,2])){
        if (is.na(the.Ct.values[sample,Detector]))
          warning(paste(sample,"is a reference sample but has no Ct value for gene",Detector))
      }
    
    
    values.of.reference.gene <- the.Ct.values[,housekeepingGene]
    if (any (is.na(values.of.reference.gene))){
      b <- names(values.of.reference.gene)[is.na(values.of.reference.gene)]
      warning(paste("There is/are no Ct values of the reference gene for the following sample/s:",paste(b,collapse=",")))}
    
    ## the delta CT value
    dCt <- the.Ct.values - values.of.reference.gene
    
    ## the delta delta CT value
    ref2a  <- dCt[calibrationSample,,drop=FALSE]
    ref2b <- apply(ref2a,2,mean,na.rm=TRUE)
    ddCt <- t (t(dCt)-ref2b)

    ## the level
    the.level <- apply(ddCt,c(1,2),function(x) 2^(-x))

    ## putting the stuff together
    cali <- rownames(ddCt)==calibrationSample
    names(cali) <-rownames(ddCt) 

    a  <- new("phenoData",
              pData=data.frame(Calibrator=cali),
              varLabels=list(Calibrator="given by user"))
    result <- new("eSet",
                  eList=list(exprs=t(the.level), Ct=t(the.Ct.values),
                    dCt=t(dCt), ddCt=t(ddCt), MAD=t(the.mad.values),
                    Difference=t(the.difference.values),
                    numberNA=t(number.of.na),
                    Plate=t(the.plate)), phenoData=a)
 
    if( !is.null(sampleInformation)) {
      if( !("Sample" %in% colnames(pData(sampleInformation))))
        stop("Your phenoData must contain a column named 'Sample'.")
      the.match <- match(rownames(pData(result)),as.character(pData(sampleInformation)$Sample))
      pData(result) <- cbind(pData(result),
                             pData(sampleInformation)[the.match, colnames(pData(sampleInformation))!="Sample"])
      phenoData(result)@varLabels<- c(varLabels(phenoData(result)),
                                      varLabels(sampleInformation)[names(varLabels(sampleInformation))!="Sample"])
    }
               
    return(result)},
                      
warning = function(x){
  ww <- file(filename,open="a+")
  writeLines(as.character(x),con=ww)
  close(ww)})
}

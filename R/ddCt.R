ddCt <- function(raw.table,
                 calibrationSample,
                 housekeepingGene,
                 sampleInformation=NULL,
                 filename="warning.output.txt",
                 type="mean"){
  withCallingHandlers({
    if (! all(c("Ct","Sample","Detector","Platename")%in% colnames(raw.table)))
      stop ("Your .txt file must include columns with the following names : 'Ct','Sample','Detector','Platename'.")
 
    aaa <- raw.table$Ct
    bbb <- raw.table$Platename
    reduced.set <- raw.table[,c("Sample","Detector")]

    if (! housekeepingGene %in% reduced.set[,2])
      stop(paste("Your reference gene",housekeepingGene,"is not part of your your .txt file."))
    if (! all(calibrationSample %in% reduced.set[,1]))
      stop("At least one of your reference samples is not part of your .txt file.")
    if (! type %in% c("median","mean"))
      stop("Type must be median or mean!")

    ## some local functions
    the.difference <- function(x) return (ifelse (any(is.na(x)), NA, max(diff(sort(x)))))
    sum.na         <- function(x) {sum(is.na(x))}
    unique.plate   <- function(x) {
      if (length(unique(x))!=1) warning(paste("g-s comb. on more than one plate:",paste(unique(x),collapse=",")))
      return (unique(x)[1])}
    sqr <- function(x) x*x
                
    number.of.na          <- tapply(aaa,reduced.set,sum.na)               # Nr. der nicht ausgewerteten Punkte
    number.of.all         <- tapply(aaa,reduced.set,length)               # alle Punkte
    if (type=="median"){
      the.Ct.values         <- tapply(aaa,reduced.set,median,na.rm=TRUE)    # der Median der Triplets
      error.Ct.mad         <- tapply(aaa,reduced.set,mad,na.rm=TRUE,con=1) # der MAD in einem Triplet
      error.Ct             <- error.Ct.mad/sqrt(number.of.all)
    } else{
      the.Ct.values         <- tapply(aaa,reduced.set,mean,na.rm=TRUE)      # der Mean der Triplets
      error.Ct.sd           <- tapply(aaa,reduced.set,sd,na.rm=TRUE)     # die SD in einem Triplet
      error.Ct              <- error.Ct.sd/sqrt(number.of.all) 
    }
    the.difference.values <- tapply(aaa,reduced.set,the.difference)       # Verhältnis lange zu kurze Strecke
    the.plate             <- tapply(bbb,reduced.set,unique.plate)
    
    for (sample in calibrationSample)
      for( Detector in unique(reduced.set[,2])){
        if (is.na(the.Ct.values[sample,Detector]))
          warning(paste(sample,"is a reference sample but has no Ct value for gene",Detector))
      }
    
    
    Ct.of.reference.gene       <- the.Ct.values[,housekeepingGene]
    error.Ct.of.reference.gene <- error.Ct[,housekeepingGene]
    if (any (is.na( Ct.of.reference.gene))){
      b <- names(Ct.of.reference.gene)[is.na(Ct.of.reference.gene)]
      warning(paste("There is/are no Ct values of the reference gene for the following sample/s:",paste(b,collapse=",")))}
    
    ## the delta CT value and errors
    dCt       <- the.Ct.values - Ct.of.reference.gene
    error.dCt <- sqrt(sqr(error.Ct) + sqr(error.Ct.of.reference.gene))
    
    
    ##erst mal fŽür nur ein Callibration Sample !!!!!
    dCt.calibration.sample       <- dCt[calibrationSample,]
    error.dCt.calibration.sample <- error.dCt[calibrationSample,]
    
    ## the delta delta CT value and errors
    ddCt <- t (t(dCt)-dCt.calibration.sample)
    error.ddCt <- t( sqrt(sqr(t(error.dCt)) + sqr(error.dCt.calibration.sample)))

    ## the delta delta CT value
    ## ref2a  <- dCt[calibrationSample,,drop=FALSE]
    ## ref2b <- apply(ref2a,2,mean,na.rm=TRUE)
    ## ddCt <- t (t(dCt)-ref2b)

    ## the levelfuntion
    levelfkt <- function(x) 2^(-x)
    
    ## level values and errors 
    the.level       <- apply(ddCt,c(1,2),levelfkt)
    the.level.error <- log(2) * the.level *  error.ddCt
    
    ## the plusfehler
    ##plusf  <-  - the.level + levelfkt(ddCt-dSD)
    ## theminusfehler
    ##minusf <- the.level - levelfkt(ddCt+dSD)

    ## putting the stuff together
 
    require(Biobase)
    cali <- rownames(ddCt)==calibrationSample
    names(cali) <-rownames(ddCt) 
    
    a  <- new("phenoData", pData=data.frame(Calibrator=cali), varLabels=list(Calibrator="given by user"))
    result <- new("eSet", eList=list (exprs=t(the.level),
                            level.err= t(the.level.error),
                            Ct=t(the.Ct.values),
                            Ct.error=t(error.Ct),
                            dCt=t(dCt),
                            dCt.error=t(error.dCt),
                            ddCt=t(ddCt),
                            ddCt.error=t(error.ddCt),
                            Difference=t(the.difference.values),
                            numberNA=t(number.of.na),
                            number  = t(number.of.all),
                            Plate=t(the.plate)),
                  phenoData=a)
    
    if (! is.null(sampleInformation)) {
      if( !("Sample" %in% colnames(pData(sampleInformation))))
        stop("Your phenoData must contain a column named 'Sample'.")
      the.match <- match(rownames(pData(result)),as.character(pData(sampleInformation)$Sample))
      pData(result) <- cbind(pData(result),
                             pData(sampleInformation)[the.match,colnames(pData(sampleInformation))!="Sample"])
      phenoData(result)@varLabels<- c(varLabels(phenoData(result)),
                                      varLabels(sampleInformation)[names(varLabels(sampleInformation))!="Sample"])
    }
               
    return(result)},
  warning = function(x){ww <- file(filename,open="a+");writeLines(as.character(x),con=ww);close(ww)})
}

miniplot <- function(result,sample)
  { e <- eList(result)
    nn <- colnames(eList(result)[[1]])
    ww <- rownames(eList(result)[[1]])
    elevel <- e$exprs[,nn==sample]
    fehler <- e$level.err[,nn==sample]
    plot(elevel, pch=16)
    #axis(2, las=2)
    #axis(1,1:length(ww),labels=ww,las=2,cex.axis=0.5)
    points(elevel + fehler , pch="-")
    points(elevel - fehler, pch="-")
  }

ddCt <- function(raw.table,calibrationSample,housekeepingGenes,type="mean",sampleInformation=NULL,filename="warning.output.txt"){
 withCallingHandlers({
 if (! all(c("Ct","Sample","Detector","Platename")%in% colnames(raw.table))) stop ("Your table must include columns with the following names : 'Ct','Sample','Detector','Platename'.")

 require(Biobase) 
 aaa <- raw.table$Ct
 bbb <- raw.table$Platename
 reduced.set <- raw.table[,c("Sample","Detector")]

 if (!all(housekeepingGenes %in% reduced.set[,2]))   stop("Not all of your housekeeping genes are in your table")
 if (! all(calibrationSample %in% reduced.set[,1]))  stop("At least one of your reference samples is not in your table.")
 if (! type %in% c("median","mean"))                 stop("Type must be median or mean!")


 the.difference <- function(x) return (ifelse (any(is.na(x)), NA, max(diff(sort(x)))))
 sum.na         <- function(x) {sum(is.na(x))}
 unique.plate   <- function(x) {
                   if (length(unique(x))!=1) warning(paste("g-s comb. on more than one plate:",paste(unique(x),collapse=",")))
                   return (unique(x)[1])}

 number.of.na           <- tapply(aaa,reduced.set,sum.na)                  # number of points with NA
 number.of.all          <- tapply(aaa,reduced.set,length)                  # number of replication
 if (type=="median"){
  the.Ct.values         <- tapply(aaa,reduced.set,median,na.rm=TRUE)       # Median
  error.Ct.mad          <- tapply(aaa,reduced.set,mad,na.rm=TRUE,con=1)    # MAD 
  error.Ct              <- error.Ct.mad/sqrt(number.of.all - number.of.na )
 } else{
  the.Ct.values         <- tapply(aaa,reduced.set,mean,na.rm=TRUE)         # Mean
  error.Ct.sd           <- tapply(aaa,reduced.set,sd,na.rm=TRUE)           # SD 
  error.Ct              <- error.Ct.sd/sqrt(number.of.all - number.of.na ) # SEM 
 }
 the.difference.values  <- tapply(aaa,reduced.set,the.difference)          # ratio long distance short distance
 the.plate              <- tapply(bbb,reduced.set,unique.plate)

 
## warning messages if a reference sample or a housekeeping gene has no values


 for (sample in calibrationSample)
  for( Detector in unique(reduced.set[,2])){
   if (is.na(the.Ct.values[sample,Detector]))
     warning(paste("No value for gene",Detector,"in ref. sample",sample))
  }

 for (Detector in housekeepingGenes)
  for( sample in unique(reduced.set[,1])){
   if (is.na(the.Ct.values[sample,Detector]))
     warning(paste("No value for housekeeping gene",Detector,"in sample",sample))
  }



# if (any (is.na( Ct.of.reference.gene))){
#   b <- names(Ct.of.reference.gene)[is.na(Ct.of.reference.gene)]
#   warning(paste("There is/are no Ct values of the reference gene for the following sample/s:",paste(b,collapse=",")))}

#########################################################
# Ct and error calculation for the housekeeping gene(s) #
#########################################################

 

 Ct.of.reference.gene       <- rowMeans(the.Ct.values[,housekeepingGenes,drop=FALSE])
 all.housekeeping.error     <- error.Ct[,housekeepingGenes,drop=FALSE]  
 not.na.per.row             <- length(housekeepingGenes) - apply(all.housekeeping.error,1,sum.na)
 error.Ct.of.reference.gene <- 1/not.na.per.row * sqrt(rowSums((all.housekeeping.error)^2))


#################################
# the delta CT value and errors #
#################################
 
 dCt         <- the.Ct.values - Ct.of.reference.gene
 
 if(length(housekeepingGenes)==1){
   error.dCt <- sqrt((error.Ct)^2 + (error.Ct.of.reference.gene)^2)
   error.dCt[ ,colnames(the.Ct.values) == housekeepingGenes] <- 0
 }else{
   the.hkg.c   <- 1 - (colnames(the.Ct.values) %in% housekeepingGenes) * 2/not.na.per.row
   red.error   <- t( t((error.Ct)^2) *the.hkg.c)
   error.dCt   <- sqrt(red.error + (error.Ct.of.reference.gene)^2)
}
#########################################################
# dCt and error calculation for the ref samples         #
#########################################################


 dCt.calibration.sample        <- colMeans(dCt[calibrationSample,,drop=FALSE])
 all.ref.sample.error          <- error.dCt[calibrationSample,,drop=FALSE]  
 not.na.per.col                <- length(calibrationSample) - apply(all.ref.sample.error,2,sum.na)
 error.dCt.calibration.sample  <- 1/not.na.per.col * sqrt(colSums((all.ref.sample.error)^2))

#######################################
# the delta delta CT value and errors #
#######################################
 
 ddCt          <- t (t(dCt)-dCt.calibration.sample)

 if(length(calibrationSample)==1){
   error.ddCt <- t( sqrt((t(error.dCt))^2 + (error.dCt.calibration.sample)^2))
   error.ddCt[rownames(the.Ct.values) == calibrationSample,] <- 0
 }else{  
   the.ref.c     <- 1 - (rownames(the.Ct.values) %in% calibrationSample) * 2/not.na.per.col
   red.error2    <- (error.dCt)^2 *the.ref.c
   error.ddCt    <- t(sqrt(t(red.error2) + (error.dCt.calibration.sample)^2))
 }

##########################
# level values and error #
##########################

 levelfkt <- function(x) 2^(-x)

 the.level       <- apply(ddCt,c(1,2),levelfkt)
 the.level.error <- log(2) * the.level *  error.ddCt

##############################
# putting the stuff together #
##############################

 cali <- rownames(ddCt) %in% calibrationSample
 names(cali) <-rownames(ddCt) 

     a  <- new("phenoData", pData=data.frame(Calibrator=cali), varLabels=list(Calibrator="given by user"))
 result <- new("eSet",      eList=list(exprs = t(the.level),
                                  level.err  = t(the.level.error),
                                  Ct         = t(the.Ct.values),
                                  Ct.error   = t(error.Ct),
                                  dCt        = t(dCt),
                                  dCt.error  = t(error.dCt),
                                  ddCt       = t(ddCt),
                                  ddCt.error = t(error.ddCt),
                                  Difference = t(the.difference.values),
                                  numberNA   = t(number.of.na),
                                  number     = t(number.of.all),
                                  Plate      = t(the.plate)),
                             phenoData  = a)
 
 if (! is.null(sampleInformation)) {
   if( !("Sample" %in% colnames(pData(sampleInformation)))) stop("Your phenoData must contain a column named 'Sample'.")
   the.match <- match(rownames(pData(result)),as.character(pData(sampleInformation)$Sample))
   pData(result) <- cbind(pData(result),pData(sampleInformation)[the.match,colnames(pData(sampleInformation))!="Sample"])
   phenoData(result)@varLabels<- c(varLabels(phenoData(result)),varLabels(sampleInformation)[names(varLabels(sampleInformation))!="Sample"])
 }


 return(result)},
  warning = function(x){ww <- file(filename,open="a+");writeLines(as.character(x),con=ww);close(ww)})
}

readSDM <- function (files,withoutPath=TRUE){
  Ctvalues <- c()
  for (file.name in  files){
    x <- scan(file=file.name,what="character",sep="\n")
    CTstart <- grep("Well",x)
    CTend   <- grep("Summary",x)
    if(length(CTstart)==0 | length(CTend)==0) stop("Your file does not seem to be a .sdm file!")
    
    number.of.skips <- CTstart
    number.of.rows  <- CTend - number.of.skips -1
    w <- read.table(file.name,sep="\t",nrows=number.of.rows,skip=number.of.skips,header=TRUE,as.is=TRUE)
    if (! all (c("Ct","Detector","Sample")%in% colnames(w))) stop("Your file does not contain the columns 'Ct','Detector' and 'Sample.")
    w <- w[,c("Sample","Detector","Ct")]
    if (withoutPath){a <- unlist(strsplit(file.name, .Platform$file.sep)); file.name <- a[length(a)]}
    Platename <- rep(file.name, nrow(w)) 
    Ctvalues <- rbind(Ctvalues,cbind(w,Platename))}
  ow <- getOption("warn")
  options(warn=-1)
  Ctvalues[,"Ct"] <- as.numeric( Ctvalues[,"Ct"])
  options(warn=ow)
  Ctvalues[,"Platename"] <- as.character( Ctvalues[,"Platename"])
  
  return(Ctvalues)}


ddCt <- function(raw.table,name.referenz.sample,name.referenz.gene){

 the.warning <- c()
 if (! all(c("Ct","Sample","Detector")%in% colnames(raw.table))) stop ("Your .txt file must include columns with the following names : 'Ct','Sample','Detector'")
  
 aaa <- as.numeric(as.character(raw.table$Ct))
 reduced.set <- raw.table[,c("Sample","Detector")]

 if (! name.referenz.gene %in% reduced.set[,2])        stop(paste("Your reference gene",name.referenz.gene,"is not part of your your .txt file."))
 if (! all(name.referenz.sample %in% reduced.set[,1])) stop("At least one of your reference samples is not part of your .txt file.")


 for (sample in name.referenz.sample)
  for( Detector in levels(reduced.set[,2])){
   if (! (sample %in% reduced.set[reduced.set[,2]==Detector,1])) the.warning <- c(the.warning,paste(sample,"is a reference sample but is not present for gene",Detector))
  }

 number.of.na          <- tapply(aaa,reduced.set,sum.na)               # Nr. der nicht ausgewerteten Punkte
 the.Ct.values         <- tapply(aaa,reduced.set,median,na.rm=TRUE)    # der Median der Triplets
 the.mad.values        <- tapply(aaa,reduced.set,mad,na.rm=TRUE,con=1) # der MAD in einem Triplet
 the.difference.values <- tapply(aaa,reduced.set,the.difference)       # Verhältnis lange zu kurze Strecke

 values.of.reference.gene <- the.Ct.values[,name.referenz.gene]
 if (any (is.na(values.of.reference.gene))){
   b <- names(values.of.reference.gene)[is.na(values.of.reference.gene)]
  the.warning <- c(the.warning,paste("There is/are no Ct values of the reference gene for the following sample/s:",paste(b,collapse=",")))}

# der delta CT Wert: neuer Mittelwert - neuer Mittelwert von einem Referenzgen
 dCt <- the.Ct.values - values.of.reference.gene

# der Delta Delta CT Wert : Delta CT - gemitteltem Delta CT Wert von Referenzsamplen
 ref2a  <- dCt[name.referenz.sample,,drop=FALSE]
 ref2b <- apply(ref2a,2,mean,na.rm=TRUE)
 ddCt <- t (t(dCt)-ref2b)

# der Level
 the.level <- apply(ddCt,c(1,2),function(x) 2^(-x))

 the.result.of.all <- cbind(convert(the.Ct.values),
                            convert(dCt)[,3],
                            convert(ddCt)[,3],
                            convert(the.level)[,3],
                            convert(the.mad.values)[,3],
                            convert(the.difference.values)[,3],
                            convert(number.of.na)[,3])

 colnames(the.result.of.all) <- c("Sample","Detector","Ct","dCt","ddCt","Level","Triplet_MAD","Long distance (for triplet)","Number_of_NA")

 return(list(all.table=the.result.of.all,levelMatrix=the.level,dCtMatrix=dCt,warn=the.warning))
}

 

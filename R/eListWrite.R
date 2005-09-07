eListWrite <- function(x,...){
   stopifnot(class(assayData(x))=="list")
   a <- data.frame(expand.grid(rownames(exprs(x)),colnames(exprs(x))))
   for (i in 1:length(assayData(x))) {b<-colnames(a)
                                  a <- cbind(a,as.vector(assayData(x)[[i]]))
                                  colnames(a) <- c(b,names(assayData(x)[i]))}
   write.table(a,...)}

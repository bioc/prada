eListWrite <- function(x,...){
   a <- data.frame(expand.grid(rownames(exprs(x)),colnames(exprs(x))))
   for (i in 1:length(eList(x))) {b<-colnames(a)
                                  a <- cbind(a,as.vector(eList(x)[[i]]))
                                  colnames(a) <- c(b,names(eList(x)[i]))}
   write.table(a,...)}

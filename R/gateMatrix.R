## ===========================================================================
## interactive locator
## ---------------------------------------------------------------------------
myLocator <- function(mymin=0,mymax=1023){
  mypoint <- locator(1)
  if (!is.null(mypoint)){
    if (mypoint$x < mymin) mypoint$x <- mymin
    if (mypoint$x > mymax) mypoint$x <- mymax
    if (mypoint$y < mymin) mypoint$y <- mymin
    if (mypoint$y > mymax) mypoint$y <- mymax
  } # end if
  return(mypoint)
}
## ===========================================================================


## ===========================================================================
## helper function to draw lines on plot
## ---------------------------------------------------------------------------
drawLine <- function(point1,point2,gatecol="red"){
  x1 <- round(point1$x)
  x2 <- round(point2$x)
  y1 <- round(point1$y)
  y2 <- round(point2$y)
  coord <- list(x=c(x1,x2),y=c(y1,y2))
  lines(coord,lty=1,lwd=2,col=gatecol)
} #drawLine
## ===========================================================================


## ===========================================================================
## wrapper for C-function 
## ---------------------------------------------------------------------------
insidePolygon <- function(data, polygon){
    # determine if a point 'mypoint' is inside a polygon given by its vertices:
    # if number of intersections between a horizontal ray emanating from
    # mypoint and edges of the polyong is odd -> mypoint inside polygon
    vertices <- matrix(as.numeric(unlist(polygon)), ncol=2, byrow=TRUE)
    data <- cbind(as.numeric(data[,1]), as.numeric(data[,2]))
    res <- .Call("inPolygon", data, vertices,
                 PACKAGE="prada")
    return(as.logical(res))
  } # insidePolygon
## ===========================================================================


## ===========================================================================
## do polygonial gating on two dimensional data
## ---------------------------------------------------------------------------
gatePoints <- function(obj,totmin=0,totmax=1023,gatecol="red",smooth=FALSE,
                       keep=TRUE, vertices=list(), comb=FALSE, add=TRUE){
  
  ##  check arguments:
  if (!is.matrix(obj)) stop("First argument is no matrix!\n")
  if (ncol(obj)!=2) stop("Data matrix does not have two columns!\n")
  if (is.null(colnames(obj))) colnames(obj) <- c("x","y")
  options(locatorBell=FALSE)
  if (is.null(totmin)) totmin <- round(min(obj))
  if (is.null(totmax)) totmax <- round(max(obj))
  mylimits <- c(totmin,totmax)

  
  # type of gate
  userAnswer <- "n"
  while(!userAnswer %in% c("p", "r"))
    userAnswer <- tolower(readline(paste("Type of gate? (P)olygon or",
                                         "(R)ectangular?")))
  rect <- ifelse(userAnswer=="r", TRUE, FALSE)

  
  # start plot
  if (smooth){
    smoothScatter(obj, xlab=colnames(obj)[1], ylab=colnames(obj)[2],
                  xlim=mylimits, ylim=mylimits)
    points(obj[prev,,drop=FALSE], col="blue", pch=1,cex=0.7)
  }else{
    plot(obj[keep,,drop=FALSE], pch=20, cex=0.5, xlab=colnames(obj)[1],
         ylab=colnames(obj)[2],
         xlim=mylimits, ylim=mylimits, col=densCols(obj)[keep])
    if(add)
      points(obj[!keep,,drop=FALSE], pch=20, cex=0.5, col="lightgray")
    if(add && length(vertices)>0){
      for(i in 1:length(vertices))
          lines(vertices[[i]], col="darkred", lwd=2)
    }
  }
  polyVertices <- c()

 
  # get input
  if(!rect){
    cat("Draw gate (left click: set point, right click: end drawing).\n")
    startpoint <- myLocator(totmin,totmax)
    polyVertices <- c(polyVertices,list(c(startpoint$x,startpoint$y)))
    temppoint <- "bla"
    lastpoint <- startpoint
    while (!is.null(temppoint)){
      temppoint <- myLocator(totmin,totmax)
      if (!is.null(temppoint)){
        polyVertices <- c(polyVertices,list(c(temppoint$x,temppoint$y)))
        drawLine(lastpoint,temppoint,gatecol)
        lastpoint <- temppoint
      }
    } # while
    drawLine(lastpoint,startpoint,gatecol)
    polyVertices <- c(polyVertices,list(c(startpoint$x,startpoint$y)))
  }else{
    cat("Draw gate (left top left and bottom right corners).\n")
    startpoint <- myLocator(totmin,totmax)
    polyVertices <- c(polyVertices,list(c(startpoint$x,startpoint$y)))
    lastpoint <- myLocator(totmin,totmax)
    p <- matrix(ncol=2, nrow=4)
    p[1,] <- unlist(startpoint)
    p[2,] <- c(lastpoint$x, startpoint$y)
    p[3,] <- unlist(lastpoint)
    p[4,] <- c(startpoint$x, lastpoint$y)
    lines(rbind(p, unlist(startpoint)), lwd=2,col=gatecol)
    polyVertices <- c(polyVertices,list(p[2,]), list(p[3,]), list(p[4,]),
                      list(p[1,]))
  }
    
  
  # process drawn box
  cat("Determining events within gate...")
  dataInBox <- insidePolygon(obj,polyVertices)
  if (smooth){
    points(obj[which(dataInBox),,drop=FALSE],pch=20,cex=0.5,col=gatecol)
  }else{
    points(obj[keep,,drop=FALSE],pch=20,cex=0.5, col=densCols(obj)[keep])
    points(obj[!(dataInBox | comb),,drop=FALSE],pch=20,cex=0.5,col="lightgray")
    lines(matrix(unlist(polyVertices), ncol=2, byrow=T), col=gatecol, lwd=2)
    if(add && length(vertices)>0){
      for(i in 1:length(vertices))
          lines(vertices[[i]], col="darkred", lwd=2)
    }
  }

  gFun <- function(x) insidePolygon(data=x, polyVertices)
  return(list(indices=dataInBox, gFun=gFun, gCol=colnames(obj),
              vertices=matrix(unlist(polyVertices), ncol=2, byrow=TRUE),
              type=ifelse(rect, "rectangle", "polygon")))
}
## ===========================================================================


## ===========================================================================
## wrapper for gatePoints to gate matrices
## ---------------------------------------------------------------------------
gateMatrix <- function(object,gate.colour="red", smooth=FALSE,
                          data.min=0,data.max=1023, name="Gate"){
  ### 0. check arguments ###
  if (!is.matrix(object))
    stop("\nfunction 'gate.matrix' can only be applied to matrices!\n")
  if (is.null(colnames(object)))
    colnames(object) <- paste("Variable",1:ncol(object),sep="")
  if(!is.character(name) || length(name)!=1)
    stop("\nmandatory parameter 'name' must be character vector of length 1") 
      
  ### 2. initialize result ###
  keepRows <- !logical(nrow(object))
  userAnswer <- "r"
  gList <- vertices <- list()
  counter <- 1
  ovars <- NULL
 
  while(userAnswer!="f"){ # if not finished
    combRows <- logical(nrow(object))
    selectedVars <- getGateVariables(object)  #prompt to select vars
    userAnswer <- "r" #make sure to enter next while loop
    while(userAnswer=="r"){ # if redo gating
      add <- (is.null(ovars) || setequal(ovars, selectedVars))
      ovars <- selectedVars
      thisGate1 <- gatePoints(object[,selectedVars,drop=FALSE],
                              totmin=data.min, totmax=data.max,
                              gatecol=gate.colour, smooth=smooth,
                              keep=keepRows, add=add)
      gate1 <- thisGate1$indices
      cat("found",sum(gate1),"\n")
      userAnswer <- "x"
      while(! userAnswer %in% c("r", "c", "p", "f"))  
        userAnswer <- tolower(readline(paste("(R)edo this gating, ",
                              "(C)ombine this gating with another one",
                              ", (P)roceed or (F)inish? ")))
    }#end while redo
    vertices <- c(vertices, list(thisGate1$vertices))
    if (sum(gate1)==0) cat("No events in drawn gate!\n") #special case no cells in gate
    else if (userAnswer %in% c("p","f")){ #if proceed or finish
      #gate is logical 'AND'
      gList[[counter]] <- new("gate", name=paste("G", counter, sep=""),
                              gateFun=thisGate1$gFun, colnames=thisGate1$gCol,
                              logic="&", type=thisGate1$type,
                              boundaries=thisGate1$vertices)               
      keepRows <- keepRows & gate1
      vertices <- list()
      counter <- counter+1
    } else if (userAnswer=="c"){ #if combination
     #gate is logical 'AND'
      logic="&"
      gList[[counter]] <- new("gate", name=paste("G", counter, sep=""),
                              gateFun=thisGate1$gFun, colnames=thisGate1$gCol,
                              logic=logic, type=thisGate1$type,
                              boundaries=thisGate1$vertices)
      combRows <-  keepRows & gate1
      userAnswer <- "a" #make sure to enter next while loop
      while (userAnswer=="a"){
        #gate is logical 'OR
        counter <- counter+1
        logic <- "|" #gate is logical 'OR'
        selectedVars <- getGateVariables(object) #prompt to select vars
        add <- setequal(selectedVars, ovars)
        userAnswer <- "r" #make sure to enter next while loop
        while (userAnswer=="r"){
          thisGate2 <- gatePoints(object[,selectedVars,drop=FALSE],
                                  totmin=data.min,totmax=data.max,
                                  gatecol=gate.colour,smooth=smooth,
                                  comb=combRows, keep=keepRows,
                                  vertices=thisGate1$vertices, add=add)
          gate2 <- thisGate2$indices 
          gList[[counter]] <- new("gate", name=paste("G", counter, sep=""),
                              gateFun=thisGate2$gFun, colnames=thisGate2$gCol,
                              logic=logic, type=thisGate2$type,
                              boundaries=thisGate2$vertices)
          cat("found",sum(gate2),"\n")
          userAnswer <- "x"
          while(! userAnswer %in% c("r", "a", "u", "f")) 
          userAnswer <- tolower(readline(paste("(R)edo this gating, (A)dd",
                                "another gating to combination, (U)se this,",
                                 "combination? ")))# prompt for input
        } #end while (userAnswer=="r")
        combRows <-  combRows | gate2
        vertices <- c(vertices, list(thisGate2$vertices))
      } #end while (userAnswer=="a")
      vertices <- list()
      keepRows <- combRows
      combRows <- TRUE
      cat("Total of",length(keepRows),"events within this gate combination.\n")
      userAnswer <- "x"
          while(! userAnswer %in% c("p", "f")) 
      userAnswer <- tolower(readline("(P)roceed or (F)inish? "))
      counter <- counter+1
    } #end else combination
  }#end while userAnswer!="f"
  #inGate[keepRows] <- TRUE
  #inPercent <- 100 * sum(inGate)/length(inGate)
  #cat("\nWithin last gates: ",inPercent,"% of all",length(inGate),"events.\n\n")
  #if(length(logic)!=1)
  #  logic <- logic[-1]
  dev.off()
  names(gList) <- sapply(gList, names)
  gset <- new("gateSet", name=name, glist=gList)
  return(gset)
}
## ===========================================================================


## ===========================================================================
## chose variables from matrix to apply gating to
## ---------------------------------------------------------------------------
getGateVariables <- function(object){
  varnames <- colnames(object)
  if(ncol(object)==2)
    return(varnames)
  selectedVar1 = selectedVar2 <- NA
  if ((.Platform$OS.type=="windows") && interactive()){
    cat("Select first variable for gating...")
    while (selectedVar1 %in% c(NA,""))
      selectedVar1 <- select.list(varnames,multiple=FALSE)
    cat(selectedVar1,"\n")
    cat("Select second variable for gating...")
    while (selectedVar2 %in% c(NA,""))
      selectedVar2 <- select.list(varnames,multiple=FALSE)
    cat(selectedVar2,"\n")
    selectedVars <- c(selectedVar1,selectedVar2)
  } else { # Unix
    cat("\nChoose two of the following variables for gating:\n")    
    print(varnames)
    while(is.na(selectedVar1)){
      varanswer1 <- readline("Variable 1 ? ")
      selectedVar1 <- varnames[as.numeric(varanswer1)]
      #if (is.na(selectedVar1))
      #  selectedVar1 <- varnames[grep(varanswer1,varnames,ignore.case=TRUE)[1]]
      } #while
    while(is.na(selectedVar2)){
      varanswer2 <- readline("Variable 2 ? ")
      selectedVar2 <- varnames[as.numeric(varanswer2)]
      #if (is.na(selectedVar2))
      #  selectedVar2 <- varnames[grep(varanswer2,varnames,ignore.case=TRUE)[1]]
    } #while
    selectedVars <- c(selectedVar1,selectedVar2)
      cat("Selected Variables:")
    print(selectedVars)
  } # else Unix
  return(selectedVars)
}
## ===========================================================================


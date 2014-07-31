gatePoints <- function(obj,totmin=0,totmax=1023,gatecol="red",smooth=FALSE){

  ##  check arguments:
  if (!is.matrix(obj)) stop("First argument is no matrix!\n")
  if (ncol(obj)!=2) stop("Data matrix does not have two columns!\n")
  if (is.null(colnames(obj))) colnames(obj) <- c("x","y")
  options(locatorBell=FALSE)

  ##  auxiliary functions:
  myLocator <- function(mymin=0,mymax=1023){
    mypoint <- locator(1)
    if (!is.null(mypoint)){
      if (mypoint$x < mymin) mypoint$x <- mymin
      if (mypoint$x > mymax) mypoint$x <- mymax
      if (mypoint$y < mymin) mypoint$y <- mymin
      if (mypoint$y > mymax) mypoint$y <- mymax
    } #if (!is.null(mypoint))
    return(mypoint)
  } #myLocator

  drawLine <- function(point1,point2,gatecol="red"){
    x1 <- round(point1$x)
    x2 <- round(point2$x)
    y1 <- round(point1$y)
    y2 <- round(point2$y)
    coord <- list(x=c(x1,x2),y=c(y1,y2))
    lines(coord,lty=1,lwd=2,col=gatecol)
  } #drawLine

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

  if (is.null(totmin)) totmin <- round(min(obj))
  if (is.null(totmax)) totmax <- round(max(obj))
  mylimits <- c(totmin,totmax)

  # start plot
  if (smooth)
    smoothScatter(obj, xlab=colnames(obj)[1], ylab=colnames(obj)[2],
                  xlim=mylimits, ylim=mylimits)
  else
    plot(obj,pch=20, cex=0.5, xlab=colnames(obj)[1], ylab=colnames(obj)[2],
         xlim=mylimits, ylim=mylimits, col=densCols(obj))

  polyVertices <- c()

  # get input
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

  # process drawn box
  cat("Determining events within gate...")
  #eturn(polyVertices)
  dataInBox <- insidePolygon(obj,polyVertices)
  if (smooth)
    points(obj[which(dataInBox),,drop=FALSE],pch=20,cex=0.5,col=gatecol)
  else
    points(obj[which(dataInBox),,drop=FALSE],pch=1,cex=0.7,col=gatecol)
  return(dataInBox)
} # gatePoints


gate.matrix <- function(object,gate.colour="red",use.smoothScatter=FALSE,
                          data.min=0,data.max=1023,max.observations=20000)
{
  ### 0. check arguments ###
  if (!is.matrix(object))
    stop("Function 'gate.matrix' can only be applied to matrices!\n")
  n.obs <- nrow(object)
  if (n.obs > max.observations){
    object <- object[1:max.observations,,drop=FALSE]
    warning(paste("gate.matrix: matrix contained more than",max.observations,
                  "observations! Only first", max.observations, "observations",
                  "were used."))
  } #  if (n.obs > max.observations)
  if (is.null(colnames(object)))
    colnames(object) <- paste("Variable",1:ncol(object),sep="")

  ### 1. define auxiliary functions ###
  getGateVariables <- function(object){
    varnames <- colnames(object)
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
        if (is.na(selectedVar1))
          selectedVar1 <- varnames[grep(varanswer1,varnames,ignore.case=TRUE)[1]]
      } #while
      while(is.na(selectedVar2)){
        varanswer2 <- readline("Variable 2 ? ")
        selectedVar2 <- varnames[as.numeric(varanswer2)]
        if (is.na(selectedVar2))
          selectedVar2 <- varnames[grep(varanswer2,varnames,ignore.case=TRUE)[1]]
      } #while
      selectedVars <- c(selectedVar1,selectedVar2)
      cat("Selected Variables:")
      print(selectedVars)
    } # else Unix
    return(selectedVars)
  } #getGateVariables

  ### 2. initialize result ###
  inGate <- rep(FALSE,nrow(object))# initialize
  keepRows <- 1:nrow(object)
  userAnswer <- "r"

  while (length(grep("^[fF]",userAnswer))==0){ # if not finished
    combRows <- numeric()
    keepobject <- object[keepRows,,drop=FALSE]
    selectedVars <- getGateVariables(object)
    userAnswer <- "r"
    while (length(grep("^[rR]",userAnswer))!=0){
      gate1 <- gatePoints(keepobject[,selectedVars,drop=FALSE],totmin=data.min,totmax=data.max,gatecol=gate.colour,smooth=use.smoothScatter)
      cat("found",sum(gate1),"\n")
      userAnswer <- readline("(R)edo this gating, (C)ombine this gating with another one, (P)roceed or (F)inish? ")
    }#while
    if (sum(gate1)==0) cat("No events in drawn gate!\n")
    else if (length(grep("^[pPfF]",userAnswer))!=0){
      keepRows <- keepRows[gate1]
      cat("Keeping only events within previous gate...\n")
    } else if (length(grep("^[cC]",userAnswer))!=0){ # combination of gates?
      combRows <-  keepRows[gate1]
      userAnswer <- "a"
      while (length(grep("^[aA]",userAnswer))!=0){
        selectedVars <- getGateVariables(object)
        userAnswer <- "r"
        while (length(grep("^[rR]",userAnswer))!=0){
          gate2 <- gatePoints(keepobject[,selectedVars,drop=FALSE],totmin=data.min,totmax=data.max,gatecol=gate.colour,smooth=use.smoothScatter)
          cat("found",sum(gate2),"\n")
          userAnswer <- readline("(R)edo this gating, (A)dd another gating to combination, (U)se this combination? ")
        } # while (length(grep("^[rR]",userAnswer))!=0)
        combRows <-  c(combRows,keepRows[gate2])
      } # while (length(grep("^[aA]",userAnswer))!=0)
      keepRows <- unique(combRows)
      cat("Total of",length(keepRows),"events within this gate combination.\n")
      userAnswer <- readline("(P)roceed or (F)inish? ")
    } # else combination of gates
  }#while
  inGate[keepRows] <- TRUE
  inPercent <- 100 * sum(inGate)/length(inGate)
  cat("\nWithin last gate: ",inPercent,"% of all",length(inGate),"events.\n\n")
  class(inGate) <- "matrixGate"
  return(inGate)
} # gate.matrix


gate.cytoFrame <- function(object,gate.colour="red",use.smoothScatter=FALSE,
                          data.min=0,data.max=1023,max.observations=20000)
{
  ### 0. check arguments ###
  if (class(object)!="cytoFrame")
    stop("Function 'gate.cytoFrame' can only be applied to objects of class 'cytoFrame'!\n")

  ### 1. call 'gate.matrix'
  result <- gate.matrix(exprs(object),gate.colour=gate.colour,
                        use.smoothScatter=use.smoothScatter,
                        data.min=data.min,data.max=data.max,
                        max.observations=max.observations)
  class(result) <- "cytoFrameGate"

  return(result)
} # gate.cytoFrame


gate.cytoSet <- function(object,gate.colour="red",use.smoothScatter=FALSE,
                          data.min=0,data.max=1023,max.observations=20000)
{
  ### 0. check arguments ###
  if (class(object)!="cytoSet")
    stop("Function 'gate.cytoSet' can only be applied to objects of class 'cytoSet'!\n")

  ### 1. prepare result ###
  n.cytoFrames <- length(object)
  result <- as.list(rep("cytoFrame not gated",n.cytoFrames))
  names(result) <- phenoData(object)$name
  if (is.null(names(result)))
    names(result) <- paste("frame",1:n.cytoFrames,sep="")
  class(result) <- "cytoSetGate"
  on.exit(return(result))

  ## CONTINUE HERE! ##
  ### 2. call 'gate.matrix'
  for (i in 1:n.cytoFrames){
    cat("Working on cytoFrame '",names(result)[i],"'...\n",sep="")
    result[[i]] <- gate.cytoFrame(object[[i]],gate.colour=gate.colour,
                                  use.smoothScatter=use.smoothScatter,
                                  data.min=data.min,data.max=data.max,
                                  max.observations=max.observations)
  } # for (i in 1:n.cytoFrames)
  #return(result)

} # gate.cytoSet


plot.matrixGate <- function(x,y,x.panels=c(1,4,5),y.panels=c(2,3,6),
                            use.smoothScatter=TRUE, limits=c(0,1023),
                            plotTitle="Gated Data",gate.colour="red",...)
{
  ### assume x to be the gate and y to be the data used for gating
  stopifnot(is.matrix(y),class(x)=="matrixGate",nrow(y)==length(x))
  threePanelPlot(y,x.panels=x.panels,y.panels=y.panels,
                 use.smoothScatter=use.smoothScatter,limits=limits,
                 addPoints=y[which(x),,drop=FALSE],addCol=gate.colour,...)
  # since x is a logical vector, we can use it to index the data matrix
  invisible(NULL)
} # plot.matrixGate


plot.cytoFrameGate <- function(x,y,x.panels=c(1,4,5),y.panels=c(2,3,6),
                               use.smoothScatter=TRUE, limits=c(0,1023),
                               plotTitle="Gated Data",gate.colour="red",
                               new.device=TRUE,verbose=TRUE,...)
{
  ### assume x to be the gate and y to be the data used for gating
  stopifnot(class(y)=="cytoFrame",class(x)=="cytoFrameGate")
  data <- exprs(y)[1:length(x),,drop=FALSE]
  # in the likely case, that observations were dropped from the cytoFrame
  #   for gating, drop them here as well

  threePanelPlot(data,x.panels=x.panels,y.panels=y.panels,
                 use.smoothScatter=use.smoothScatter,limits=limits,
                 plotTitle=plotTitle,
                 addPoints=data[which(x),,drop=FALSE],addCol=gate.colour,
                 new.device=new.device,verbose=verbose,...)
  # since x is a logical vector, we can use it to index the data matrix
  invisible(NULL)
} # plot.cytoFrameGate


plot.cytoSetGate <- function(x,y,x.panels=c(1,4,5),y.panels=c(2,3,6),
                               use.smoothScatter=TRUE, limits=c(0,1023),
                               plotTitle="Gated Data",gate.colour="red",
                               verbose=TRUE,...)
{
  ### assume x to be the gate and y to be the data used for gating
  stopifnot(class(y)=="cytoSet",class(x)=="cytoSetGate")

  n.cytoFrames <- length(y)
  frame.names  <- phenoData(y)$name
  if (is.null(frame.names))
    frame.names <- paste("frame",1:n.cytoFrames,sep="")

  ### 1. prepare plot
  par(mfrow=c(n.cytoFrames,3))

  ### 2. call 'plot.cytoFrameGate'
  for (i in 1:n.cytoFrames){
    cat("Working on cytoFrame '",frame.names[i],"'...\n",sep="")

    thisFrame <- y[[i]]
    thisFrameGate <- x[[i]]
    class(thisFrameGate) <- "cytoFrameGate"
    thisPlotTitle <- paste("Gated Frame",abbreviate(frame.names[i],15))

    plot.cytoFrameGate(thisFrameGate,thisFrame,
                       x.panels=x.panels,y.panels=y.panels,
                       gate.colour=gate.colour, plotTitle=thisPlotTitle,
                       use.smoothScatter=use.smoothScatter,
                       limits = limits, new.device=FALSE, verbose=verbose,...)

  } # for (i in 1:n.cytoFrames)

  invisible(NULL)
} # plot.cytoFrameGate



threePanelPlot <- function(data,x.panels=c(1,4,5),y.panels=c(2,3,6),
                           tot.width=15,tot.height=5.4,maxcells=20000,
                           limits=c(0,1023),remove.extremes=TRUE,
                           plotTitle="Three-Panel Plot",use.smoothScatter=TRUE,
                           palette=colorRampPalette(brewer.pal(9,"Blues")),
                           new.device=TRUE,verbose=TRUE,
                           addPoints=NULL,addCol="red",...)
{
  # check arguments
  if (verbose) cat("Checking arguments...\n")
  stopifnot(all(c(x.panels,y.panels) %in% 1:ncol(data)))

  # prepare data:
  if (verbose) cat("Preparing data...\n")
  if (remove.extremes){
    data.rowmins <- apply(data,1,min)
    data.rowmaxs <- apply(data,1,max)
    data <- data[(data.rowmins>limits[1])&(data.rowmaxs<limits[2]),,drop=FALSE]
  } # if (remove.extremes)
  ncells <- nrow(data)
  sampled.events <- sample(1:ncells,min(ncells,maxcells),replace=FALSE)
  data <- data[sampled.events,,drop=FALSE]

  # prepare plot
  if (verbose) cat("Plotting...\n")
  if (new.device){
    do.call(getOption("device"),list(width=tot.width,height=tot.height))
    par(mfrow=c(1,3),font.lab=2,mar=c(5,4,5,4))
  } # if (new.device)
  for (i in 1:3){
    if (use.smoothScatter)
      smoothScatter(data[,c(x.panels[i],y.panels[i],drop=FALSE)],
                    colramp=palette,xlim=limits,ylim=limits,...)
    else
      plot(data[,c(x.panels[i],y.panels[i],drop=FALSE)],xlim=limits,
           ylim=limits,...)
    if (all(i==2,all(!is.na(plotTitle)),!is.null(plotTitle)))
      title(main=plotTitle,cex.main=1.5,font.main=2)
    if (!is.null(addPoints))
      points(addPoints[,c(x.panels[i],y.panels[i],drop=FALSE)],col=addCol,
             ...)
  } # for (i in 1:3)

} # threePanelPlot



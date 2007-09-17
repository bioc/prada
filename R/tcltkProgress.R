progress <- function(title="processing task...", message="",
                     sub=""){

  ## check for tcltk installation and exit gracefully if none
  if(!capabilities("tcltk")){
    warning("Need tcltk for the status bar")
    return(NULL)
  }
  
  ## check arguments
  options(warn=-1)
  require(tcltk)
  if(!is.character(title) || !is.character(message)|| !is.character(sub))
    stop("All arguments must be character vectors of length 1")

  ## set up necessary variables in global environment 
  assign(".tkprogress.canceled", FALSE, .GlobalEnv)
  assign(".tkprogress.window", tktoplevel(), .GlobalEnv)
  assign(".tkprogress.iterator", tclVar("0"), .GlobalEnv)
  assign(".tkprogress.fallbackIterator", tclVar("0 %"), .GlobalEnv)
  assign(".tkprogress.labelText", tclVar("0 %"), .GlobalEnv)

  ## the tcl/tk toplevel window
  tkwm.geometry(.tkprogress.window, "280x140")
  tkwm.title(.tkprogress.window, title)
  tkconfigure(.tkprogress.window, cursor="watch")
  frame <- tkframe(.tkprogress.window, relief="groove",borderwidth=2)
  font <- tkfont.create(family="helvetica",size=13,weight="bold")
  font2 <- tkfont.create(family="helvetica",size=11,weight="bold")
  font3 <- tkfont.create(family="helvetica",size=11)

  ## the textual status indicator
  if(sub!="")
    tclvalue(.tkprogress.labelText) <<- sub
  label <- tklabel(frame, text=tclvalue(.tkprogress.labelText), font=font3)
  tkconfigure(label,textvariable=.tkprogress.labelText)

  ## the header message
  if(message!=""){
    tkgrid(tklabel(.tkprogress.window, text=message, font=font2))
    tkgrid(tklabel(.tkprogress.window, text=""))
  }

  ## use BWidgetsprogress bar widget or fallback to textual indicator
  if(as.character(tclRequire("BWidget"))!="FALSE"){
    progBar <- tkwidget(.tkprogress.window, "ProgressBar",
                        variable=.tkprogress.iterator)
    tkgrid(progBar)
    tkconfigure(frame, borderwidth=0)
    tkgrid(label)
  }else{
    itLabel <- tklabel(frame, text=tclvalue(.tkprogress.fallbackIterator), 
                       font=font)
    tkconfigure(itLabel,textvariable=.tkprogress.fallbackIterator)
    tkgrid(itLabel)
    if(sub!="")
      tkgrid(label)
  }
  tkgrid(frame)
  tkbind(.tkprogress.window, "<Destroy>", function(){
    .tkprogress.canceled <<- TRUE; killProgress()})
  options(warn=0)
}

updateProgress <- function(percentage, autoKill=FALSE, sub=""){
  ## check for tcltk installation and exit gracefully if none
  if(!capabilities("tcltk"))
    return(NULL)
  percentage <- as.integer(percentage)
  tclvalue(.tkprogress.iterator) <<- percentage
  tclvalue(.tkprogress.fallbackIterator) <<- paste(percentage, "%")
  if(sub==""){
    tclvalue(.tkprogress.labelText) <<- paste(percentage, "%")
  }else{
    tclvalue(.tkprogress.labelText) <<- sub
  }
  tcl("update", "idletasks")
  if(autoKill)
    if(percentage>=100)
      killProgress()
}

killProgress <- function(){
   ## check for tcltk installation and exit gracefully if none
  if(!capabilities("tcltk"))
    return(NULL)
  if(!.tkprogress.canceled)
    tkdestroy(.tkprogress.window)
}


#progress(message="reading in multiple files...", sub="(1 of 50 files)")
#for(i in 1:50) {rnorm(1e+05); updateProgress(i*2, autoKill=TRUE, 
#                                  sub=paste("(", i, " of 50 files)", sep=""))}

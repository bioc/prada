
progress <- function(title="processing task...", message=""){
  options(warn=-1)
  require(tcltk)
  if(!is.character(title) || !is.character(message))
    stop("All arguments must be character vectors of length 1")
  assign(".tkprogress.canceled", FALSE, .GlobalEnv)
  assign(".tkprogress.window", tktoplevel(), .GlobalEnv)
  tkwm.geometry(.tkprogress.window, "250x140")
  tkwm.title(.tkprogress.window, title)
  tkconfigure(.tkprogress.window, cursor="watch")
  assign(".tkprogress.labelText", tclVar("0"), .GlobalEnv)
  frame <- tkframe(.tkprogress.window, relief="groove",borderwidth=2)
  font <- tkfont.create(family="sansserif",size=14,weight="bold")
  font2 <- tkfont.create(family="sansserif",size=12,weight="bold")
  label <- tklabel(frame, text=tclvalue(.tkprogress.labelText), font=font)
  tkconfigure(label,textvariable=.tkprogress.labelText)
  tkgrid(tklabel(.tkprogress.window, text=""))
  if(message!=""){
    tkgrid(tklabel(.tkprogress.window, text=message, font=font2))
    tkgrid(tklabel(.tkprogress.window, text=""))
  }
  if(as.character(tclRequire("BWidget"))!="FALSE"){
    ## tcltk progress bar from BWidgets library##
    progBar <- tkwidget(.tkprogress.window, "ProgressBar",
                        variable=.tkprogress.labelText)
    tkgrid(progBar)
    tkgrid(tklabel(.tkprogress.window, text=""))
    tkconfigure(frame, borderwidth=0)
  }
  tkgrid(label, tklabel(frame, text="% done", font=font))
  tkgrid(frame)
  tkgrid(tklabel(.tkprogress.window, text=""))
  tkbind(.tkprogress.window, "<Destroy>", function(){
    .tkprogress.canceled <<- TRUE; killProgress()})
  options(warn=0)
}

updateProgress <- function(percentage, autoKill=FALSE){
  percentage <- as.integer(percentage)
  tclvalue(.tkprogress.labelText) <<- percentage
  tcl("update", "idletasks")
  if(autoKill)
    if(percentage>=100)
      killProgress()
}

killProgress <- function(){
  if(!.tkprogress.canceled)
    tkdestroy(.tkprogress.window)
}


#for(i in 1:50) {rnorm(1e+05); updateProgress(i*2, auto=T)}

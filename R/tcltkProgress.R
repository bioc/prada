
progress <- function(title="processing task...", message=""){
  require(tcltk)
  if(!is.character(title) || !is.character(message))
    stop("All arguments must be character vectors of length 1")
  assign(".tkprogress.canceled", FALSE, .GlobalEnv)
  assign(".tkprogress.window", tktoplevel(), .GlobalEnv)
  tkwm.geometry(.tkprogress.window, "250x120")
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
  tkgrid(label, tklabel(frame, text="% done", font=font))
  tkgrid(frame)
  tkgrid(tklabel(.tkprogress.window, text=""))
  tkbind(.tkprogress.window, "<Destroy>", function(){ .tkprogress.canceled <<- TRUE; killProgress()})
}

updateProgress <- function(percentage, autoKill=FALSE){
  percentage <- as.integer(percentage)
  tclvalue(.tkprogress.labelText) <<- percentage
  if(autoKill)
    if(percentage>=100)
      killProgress()
}

killProgress <- function(){
  if(!.tkprogress.canceled)
    tkdestroy(.tkprogress.window)
}


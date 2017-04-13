check.Wj <-
function(Wj,doStop=TRUE){
  ## Checks to see if Wj has proper names.
  isBad <- FALSE
  if(is.null(names(Wj))){
    isBad <- TRUE
  } else {
    if(sum(names(Wj)!=paste("C",1:length(Wj)+1,sep=""))!=0)
      isBad <- TRUE
  }
  if(isBad & doStop){
    stop("Wj should be have names C2,C3, etc.")
  } else
    return(!isBad)
}

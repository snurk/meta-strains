dPhidW <-
function(Xj,Rj,Qu,Wj,doDebug=FALSE){
  ## will compute derivative of Phi with respect to W
  ## by chain rule: sum_{c'} (dPhidPcp*dPdW)
  result <- list()
  checked <- check.Wj(Wj)
  pj <- W2P(Wj) ##Pj without the first entry associated with normal sample.
  C <- length(pj)+1
  if(length(pj)!=length(Wj))##QC,
    stop("Length of Wj should be less than length of Pj!")
  if(C<2)
    stop("Too small C (clone number)!")
  ##
  gradient <- c()
  d1 <- c()
  d2 <- list()
  ## for improvement in efficiency by moving d1 & d2 out of loops
  for(cp in 2:C){
    d1[cp] <- dPhidPcp(Xj=Xj,Rj=Rj,Qu=Qu,Pj=c(0.5,pj),cp=cp)
    ##^ The first entry of Pj doesn't contribute.
  }
  for(cee in 2:C){
    d2[[cee]] <- dPdW(Wj=Wj,cee=cee)
  }
  for(cee in 2:C){
    ceeChar <- as.character(cee)
    gradient[ceeChar] <- 0
    for(cp in 2:C){
      ##^ -cp: c', the index for Z_{i,c'}. See the paper.
      gradient[ceeChar] <- gradient[ceeChar] +
        d1[cp]*(d2[[cee]][paste("C",cp,sep="")])
      if(doDebug)
        browser()
    }
  }
  ##
  result[["pj"]] <- pj ## of length 1 LESS than Pj.
  ##QC:
  if(length(gradient)!=length(Wj))
    stop("The length of gradient should be the same as input Wj!")
  result[["dPhidWOut"]] <- gradient
  return(result)
}

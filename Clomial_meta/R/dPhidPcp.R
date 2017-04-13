dPhidPcp <-
function(Xj,Rj,Qu,Pj,cp,doDebug=FALSE){
  ## will compute partial derivative of Phi with respect to P_{c',j}.
  ## -Pj is the j'th column of P with length C.
  ## -Muc is the c'th column of Mu.
  ## -cp: c', the index for Z_{i,c'}. See supplement.pdf.
  C <- length(Pj)
  derivative <- 0
  for( state in 2:(2^C) ){ ## for all possible values of a row of Z, expect 0s.
    zee <- decimal2binary(state,C) ## binary representation of state-1
    if(zee[cp]==0)
      next ## 0*x=0
    #piBinom <- (zee/2)%*%Pj ## /2 because of heterozygosity,
    ## Note: maybe derivative is twice more or less, but in the same direction.
    piBinom <- (zee)%*%Pj
    
    
    piBinom[piBinom>1] <- 1
    XOverpiFraction <- Xj/piBinom
    RXOverpiFraction <- (Rj-Xj)/(1-piBinom)
    fraction <- XOverpiFraction-RXOverpiFraction
    derivBinom <- Qu[,state]*fraction
    summand <- sum(derivBinom) ## sum over i 
    ##
    if(!is.nan(summand))
      derivative <- derivative + 1*summand ## zee[cp]=1
    if(doDebug)
      browser()
  }
  #derivative <- derivative /2  ## /2 because of heterozygosity
  result <- derivative 
  return(result)
}

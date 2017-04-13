Phi <-
function(Xj,Rj,Qu,Pj,doDebug=FALSE){
  ## will return Expected_Z[log(Pr(X|Z,Pj))]
  ## -Rj will be Dc[,j].
  ## -Xj will be Dt[,j].
  result <- 0
  C <- length(Pj)
  if(ncol(Qu)!=2^C)
    stop("Not compatible inputs!")
  for( state in 2:(2^C) ){ ## for all possible values of a row of Z,
    ## state=1 => zee=rep(0) => Pr=0 => safe to ignore
    zee <- decimal2binary(state,C) ## binary representation of state-1
    
    #piBinom <- (zee/2)%*%Pj ## /2 because of heterozygosity
    piBinom <- (zee)%*%Pj
    
    piBinom[piBinom>1] <- 1
    summand <- dbinom(Xj,Rj,piBinom,log=TRUE) ## it's a vector
    s1 <- summand
    if(doDebug)
      browser()
    ##summand <- summand+ Xj*log(piBinom)
    ##summand <- summand+ (Rj-Xj)*log(1-piBinom)
    summand <- Qu[,state]*summand
    summand[Qu[,state]==0] <- 0 ## NaN happens when Qu=piBinom=0.
    summand <- sum(summand)
    result <- result + summand
  }
  return(result)
}

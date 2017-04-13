Clomial.iterate <-
function(Dt, Dc, Mu, P, maxIt=100, U=NULL, PTrue=NULL, 
                            llCutoff=10^(-3),
                            computePFunction=compute.P.reparam,
                            doSilentOptim=TRUE,
                            doTalk=TRUE,doLog=TRUE,debug=FALSE,
                            noiseReductionRate=0.01,
                            fliProb=0.05,
                            conservative=TRUE){
  ## TODO: Add a flag such that if U and PTrue are present, 
  ## then Clomial.reorder() function (to be written)
  ## will be called and the error for each iteration
  ## is computed. 
  ## Iterates over EM.
  ## Returns a record of Mus and Ps for every iteration, the final Mu and P,
  ## and the likelihood after each iteration
  ## maxIt: maximum number of iterations allowed before returning the results
  ## Dc: matrix of coverage
  ## Dt: matrix of mutated allele counts
  ## Mu, P: initialized parameter matrices to start iterating on
  ## U, PTrue: True values for matrices. They're not even used by this function.
  ##^ Used for debugging purposes. 
  ## noiseReductionRate: The value which will be deduced from the noise
  ## in each iteration
  ## fliProb: The initial value for the noise which will affect Mu by
  ## flipping some bits randomly, i.e., 1-Mu_{i,j}.
  ##^ Set to 0 to disable  
  ## conservative: Injects the noise only if the likelihood improves after
  ## an iteration.
  if(maxIt<3){
    stop("It is recommended to have at least 3 iterations (maxIt>2) !")
  }
  N <- nrow(Mu)
  C <- ncol(Mu)
  S <- ncol(P)
  startTime <- Sys.time()
  if(is.infinite(noiseReductionRate)){
    fliProb <- 0
  }
  ##
  ## d2b is a list indexed by a decimal value whose element is
  ## the binary representation of that decimal value with length C.
  d2b <- list()
  for(i in 1:2^C){
    d2b[[i]] <- decimal2binary(i, C)
  }
  ##
  ## A state of z is a binary vector of length C that could be a row of Z.
  ## ustates[c] is a list of size 2^(C-1) of decimal representations of 
  ## the states of z where zc = 1. There are C such lists in ustates.
  ustates <- compute.unity.states(C)
  ##
  k <- 1 ## number of iterations
  log.likelihoods <- c(-Inf)
  mus <- list()
  mus[["1"]] <- Mu
  ps <- list()
  ps[["1"]] <- P
  qs <- list()
  while(k < maxIt){
    if(doTalk){
      print(paste("Iteration:",k))
      print(paste("Flipping probability:",fliProb))
    }
    ## Add noise?
    MuNoisy <- Mu
    if(fliProb > 0){
      flips <- sample(x=0:1,size=length(Mu),replace=TRUE,
                      prob=c(1-fliProb,fliProb))
      MuNoisy <- (flips)*(1-Mu)+(1-flips)*Mu
      ## No row of all 0s,
      rows0 <- rowSums(MuNoisy[,-1])==0
      MuNoisy[rows0,] <- Mu[rows0,] ## Don't change these rows
      MuNoisy[,1] <- Mu[,1] ## Don't change the first column
      ## Update noise:
      fliProb <- max((fliProb-noiseReductionRate),0)
    }
    ##
    computedQ <- compute.q(Dc,Dt,Mu,P,doLog=doLog)
    computedQNoisy <- compute.q(Dc,Dt,MuNoisy,P,doLog=doLog)
    if(debug){
      print("P");print(P)
      print("Mu");print(Mu)
      print("computedQ:");print(computedQ)
    }
    Mu <- compute.mu(computedQ$normalizedUnLog,ustates,N,C)
    MuNoisy <- compute.mu(computedQNoisy$normalizedUnLog,ustates,N,C)
    if(debug)
      print(Mu)
    ##Q=exp(qu)
    computedQ <- compute.q(Dc,Dt,Mu,P,doLog=doLog)
    computedQNoisy <- compute.q(Dc,Dt,MuNoisy,P,doLog=doLog)
    ##^ extra qu update.
    P <- computePFunction(C=C,S=S,N=N,Dt=Dt,Dc=Dc,ustates=ustates,
                          qu=computedQ$normalizedUnLog,
                          d2b=d2b,P=P,Mu=Mu,PTrue=PTrue,U=U,debug=debug,
                          doSilentOptim=doSilentOptim,doTalk=doTalk)$optimumP
    PNoisy <- computePFunction(C=C,S=S,N=N,Dt=Dt,Dc=Dc,ustates=ustates,
                               qu=computedQNoisy$normalizedUnLog,
                               d2b=d2b,P=P,Mu=MuNoisy,PTrue=PTrue,U=U,
                               debug=debug,
                               doSilentOptim=doSilentOptim,
                               doTalk=doTalk)$optimumP
    ##
    ll <- Clomial.likelihood(Dc, Dt, Mu, P)$ll
    llNoisy <- Clomial.likelihood(Dc, Dt, MuNoisy, PNoisy)$ll
    ## Add noise only if is beneficial.
    if( llNoisy>(ll+llCutoff) | !conservative){
      if(doTalk)
        print(paste("Noise was applied; from:",ll,"to:",llNoisy))
      Mu <- MuNoisy
      P <- PNoisy
      ll <- llNoisy
    }
    log.likelihoods[k] <- ll
    if(debug){
      print(k); print(log.likelihoods[k])
    }
    ##
    ## Check for convergence		
    if( k > 1 ){
      ## Do at least two iterations before stopping
      likelihood.ratio <- abs((log.likelihoods[k] - log.likelihoods[k-1])
                              /log.likelihoods[k-1])
      if(doTalk)
        print(paste("LL ratio:",likelihood.ratio))
      ##
      if( log.likelihoods[k-1] != -Inf & likelihood.ratio < llCutoff ){
        break
      } 
    }
    k <- k+1
    if(debug)
      print(log.likelihoods)
    ##    
  }##End while(k < maxIt).
  ##
  endTime <- Sys.time()
  output <- list()
  output[["Qs"]] <- qs  #### history of qs
  output[["Ps"]] <- ps  #### history of ps
  output[["Mus"]] <- mus #### history of Mus
  output[["Mu"]] <- Mu ####optimized Mu (not reordered)
  output[["P"]] <- P #### optmized P (not reordered)
  output[["llCutoff"]] <- llCutoff
  output[["maxIt"]] <- maxIt
  output[["LRatio"]] <- likelihood.ratio
  output[["Likelihoods"]] <- log.likelihoods
  output[["fliProb"]] <- fliProb
  output[["endTime"]] <- endTime
  output[["timeTaken"]] <- endTime-startTime
  return(output)
}

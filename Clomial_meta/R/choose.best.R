choose.best <-
function(models,U=NULL,PTrue=NULL,compareTo=NULL,upto="All", doTalk=FALSE){
  ## Chooses the best iteration out of the set
  ## -compareTo: The try with which each try is compared.
  ## are considered and P will be updated accordingly.
  ## The pruned model with the best likelihood
  ## will be returned. I need to write prune.genotype() function to do this.
  ## which make phyl tree
  ## and returned the one with the best likelihood.
  ## -models can be a list of outputs produced by iterate(), 
	## or it can be a vector of characters
  ## each of which pointing to a file containing the output for a random seed.
  ## 
  maxLikelihood <- -Inf
  bestIndex <- 1
  bestModel <- "none"
  model <- NULL ## Will be loaded.
  ## The error when truth is known.
  err <- list();err[["U"]] <- c();err[["P"]] <- c()
  ## The difference between a try and the best one.
  diffBest <- list(); diffBest[["U"]] <- c(); diffBest[["P"]] <- c()
  diffBest[["discretizedU"]] <- c();
  Li <- c()
  seconds <- c()
  if(upto=="All"){
    inds <- 1:length(models)
  }else{
    inds <- 1:upto
  }
  if(doTalk)
    print(paste("Choosing the best out of",length(inds),"..."))
  for( JInd in inds){
    J <- as.character(JInd)
    if(JInd %% 1000==0 & doTalk)
      print(JInd)
    Li[J] <- NA
    seconds[J] <- NA
    err[["U"]][J] <- NA
    err[["P"]][J] <- NA
    if(class(models)=="list"){
      outI <- models[[JInd]]
    }else{ ## models is the vector of saved results.
      binomTries <- NULL
      load(models[JInd]) ## binomTries is loaded.
      outI <- binomTries[[1]]
    }
    if( class(outI) != "try-error" ){
      LLI <- tail( outI[["Likelihoods"]], n = 1 )
      if( LLI > maxLikelihood ){
        bestIndex <- JInd
        maxLikelihood <- LLI
        bestModel <- outI
      }##End if.
      Li[J] <- LLI
      seconds[J] <- as.numeric(outI$timeTaken, unit="secs") 
      if(!is.null(U) & !is.null(PTrue)){
        computedError <- compute.errors(Mu=outI$Mu,U=U, P=outI$P, PTrue=PTrue)
        err[["U"]][J] <- computedError$UError
        err[["P"]][J] <- computedError$PErrorAbsolute
      }
      ##
      if(!is.null(compareTo)){
        ## How different is each try from the best?
        if(class(models)=="list"){
          compareToTry <- models[[compareTo]]
        } else {
          binomTries <- NULL
          load(models[compareTo]) ## binomTries is loaded.
          compareToTry <- binomTries[[1]]          
        }
        computedDiff <- compute.errors(Mu=outI$Mu,U=compareToTry$Mu,P=outI$P,
                                       PTrue=compareToTry$P)
        diffBest[["U"]][J] <- computedDiff$UError
        diffBest[["discretizedU"]][J] <- computedDiff$discretizedUError
        diffBest[["P"]][J] <- computedDiff$PErrorAbsolute
      }
    }##End if.    
  }##End for.
  result <- list(err=err,Li=Li,bestInd=bestIndex,bestLi=Li[bestIndex],
                 comparison=diffBest,bestModel=bestModel,seconds=seconds)
  return(result)
}

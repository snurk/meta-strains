compute.errors <-
function(Mu,U,P,PTrue){
  ## Finds the permutation of U that best matches Mu and returns the errors of
  ## Mu and P corresponding to that permutation
  ## INPUTS,
  ## P and Mu: inferred.
  ## PTrue and U: the truth.
  ##
  C <- ncol(Mu)
  if(C>7){
    stop("Too many clones to consider all permutations!")
  }
  allPermS <- rbind( 1:C, permute::allPerms(C))
  reordered <- reorder.clones(U=U,PTrue=PTrue,Mu=Mu,P=P, allPermS=allPermS)
  Mu <- reordered[["Mu"]]
  P <- reordered[["P"]]
  ##  
  UError <- sum(abs(Mu-U)) /length(Mu)
  discUError <- sum(abs(round(Mu)-round(U))) /length(Mu)
  PErrorAbsolute <- mean(abs(P-PTrue))
  PErrorRelative <- mean(abs(P-PTrue)/PTrue)
  ##
  errors <- list()
  errors[["UError"]] <- UError
  errors[["discretizedUError"]] <- discUError
  errors[["PErrorAbsolute"]] <- PErrorAbsolute
  errors[["PErrorRelative"]] <- PErrorRelative
  ##
  return(errors)
}

compute.mu <-
function( qu, ustates, N, C ){
  ## Returns the updated NxC mu matrix
  ## Each (i,c) element is the sum of all qu[i,'z'] such that zc == 1.
  ##
  mu <- matrix(0,N,C)
  rownames(mu) <- rownames(qu)
  for(c in 1:C){
    for(i in 1:N){
      for(z in 1:length( ustates[[c]] )){
        mu[i,c] <- mu[i,c] + qu[i, ustates[[c]][[z]] ]
      }		
    }		
  }
  mu[,1] <- matrix(0,N,1)
  colnames(mu) <- paste("C",1:C,sep="")
  return(mu)
}

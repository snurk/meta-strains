best.match <-
function(m1, m2, perms){
  ## Matches the columns of row m2 to the columns of m1. 
  ## Values should be all positive or all negative.
  ## Returns a list containing 
  ## [1] the ordering of the cols of m2 that corresponds to m1
  ## [2] the absolute-valued difference between the matched m1 and m2
  ## Brute force approach
  products <- t(m1)%*%(m2) 
  maxPerm <- -1
  maxSum <- 0
  for(p in 1:nrow(perms)){
    s <- 0
    for(i in 1:nrow(products)){
      s <- s+products[i,perms[p,i]]
    }
    if(s > maxSum){
      maxSum <- s
      maxPerm <- p
    }
  }
  if(maxPerm == -1){
    warning("best_match failed!")
  }
  return(list(perms[maxPerm,], sum(abs(m1 - m2[,perms[maxPerm,]]))))
}

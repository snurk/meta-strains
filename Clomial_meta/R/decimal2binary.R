decimal2binary <-
function(n, digits){
    ## Returns a 1 by C matrix that is the binary representation of n-1.
    ## For n> 2^C, we are looking at the remainder of n by 2^C.
    ## Example: (n=10, C=3) =>  10 -> 2 -> 1 -> 001
    ## Maybe could be optimized but not needed for our code.
  n <- n-1
  answer <- matrix(0,1,digits)
  for( i in 1:digits ){
    if( n %% (2^i) != 0 ){
      answer[i] <- 1
      n <- n-(2^(i-1))
    }
  }
  return( rev(answer) )
}

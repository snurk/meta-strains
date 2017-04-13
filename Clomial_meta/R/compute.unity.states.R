compute.unity.states <-
function(C){
  ## Returns a list of lists of the binary states of z
  ## where zc = 1 for each c in C.
  ## answer[c] will be a list of size 2^(C-1) of
  ## decimal representations of the states of z.
  ## This simplifies the computation of new Mus 
  answer <- list()
  for(i in 1:C){
    answer[[i]] <- list()	
  }
  for(state in 1:(2^C)){
    z <- decimal2binary(state,C)
    for(j in 1:C){
      if(z[j] == 1){
        answer[[j]][[ length(answer[[j]])+1 ]] <- state			
      }
    }
  }
  return(answer)
}

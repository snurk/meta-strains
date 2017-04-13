row.z.likelihood <-
function(Zi, mu_i){
  ## Product of P(Zi[c]) = mu_i[c] for each Zi[c]
  ## assuming a Bernoulli distribution
  C <- length(Zi)
  answer <- 1
  for(c in 1:C){
    if(Zi[c] == 1){
      answer <- answer * mu_i[c]		
    }	
    else{
      answer <- answer * ( 1-mu_i[c] )		
    }
  }
  return(answer)
}

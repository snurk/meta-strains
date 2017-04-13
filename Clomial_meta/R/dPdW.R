dPdW <-
function(Wj,cee){
  ## will compute derivative of P with respect to W,
  ## i.e, (1+sum(W)-W)/(1+sum(W))^2
  ## -Wj is a vector (probably of length C-1).
  ## -cpIsc: determines if c'=c or not.
  ## OUTPUT: the vector of gradient.
  ## QC:
  checked <- check.Wj(Wj)
  sigma <- sum(Wj)
  cpIsc <- 0*Wj
  cpIsc[paste("C",cee,sep="")] <- 1
  result <- ((1+sigma)*cpIsc-Wj) / (1+sigma)^2
  return(result)
}

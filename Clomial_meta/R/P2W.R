P2W <-
function(Pj){
  ## Computes Vj(Pj)~ Pj/(1-sum(Pj)).
  ## No log due to computational problems such as log(900)=Inf.
  Pj1 <- Pj[2:length(Pj)] ## ignore the first row of P.
  Wj <- Pj1/(1-sum(Pj1))
  return(Wj)
}

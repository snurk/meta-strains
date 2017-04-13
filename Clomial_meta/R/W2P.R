W2P <-
function(W){
    ## Computes: Pj(Wj) ~ W/(1+sum(W)).
    ## -W and P are vectors.
    P <- W/(1+sum(W))
    ##print(P)
    return(P)
  }

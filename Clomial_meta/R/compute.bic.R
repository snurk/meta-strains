compute.bic <-
function( Dc, Dt, Mu, P){
  ## Computes the BIC using:
  ## BIC = NC log(obsNum) - 2 log(L)
  ## Where n = Dc_[i,t], k = Dt_[i,t], p = and Mu_i *P[,t] 
  ## L = log-likelihood
  result <- list()
  ll <- Clomial.likelihood(Dc=Dc, Dt=Dt, Mu=Mu, P=P)$ll
  result[["ll"]] <- ll
  N <- nrow(Mu)
  C <- ncol(Mu)
  M <- ncol(P)
  obsNum <- sum(Dc)
  obsNumAverage <- obsNum/(N*M)
  aic <- -2*ll + (N*C+M*(C-1))*2
  bic <- -2*ll + (N*C+M*(C-1))*log(obsNum)
  ##
  result[["bic"]] <- bic
  result[["aic"]] <- aic
  result[["obsNum"]] <- obsNum
  return(result)
}

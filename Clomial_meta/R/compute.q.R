compute.q <-
function(Dc,Dt,Mu,P,doEnforceNormal=TRUE,doLog=TRUE){
  ## FIXME, optimization: pass d2b and avoid calling decimal2binary function
  ## Returns an Nx(2^C) matrix
  ## where each element is the Pr(Di,z|P,mu) for a particular i and z.
  ## Pr(Di,z|P,mu) = Pr(Di|z,P) * Pr(z|P,mu)
  ## Each column represents the binary string, translated into decimal,
  ## that represents a particular state of zi.
  result <- list()
  N <- nrow(Dc)
  C <- nrow(P)
  rowDataLikelihoods <- matrix(0,N,2^C)
  rowDataLogLikelihoods <- matrix(0,N,2^C)
  row.z.likelihoods <- matrix(0,N,2^C)
  qu <- matrix(0,N,2^C)
  rownames(qu) <- rownames(Dc)
  if(sum(Mu[,1])!=0 & doEnforceNormal)
    stop("The first column of Mu, corresponding to normal, is not all 0s.")
  ##
  for(i in 1:N){
    for(state in 1:(2^C)){
      z <- decimal2binary(state,C)
      ## p <- z%*%P
      ## p is \pi in the paper.
      
      #p <- (z/2)%*%P ## /2 because of heterozygosity
      p <- (z)%*%P
      
      p[p>1] <- 1
      ## Pr(X.i|Z.i,P,mu) given by product of binomials
      ##
      rowDataLikelihoods[i,state] <- prod(dbinom(Dt[i,],Dc[i,],p))
      rowDataLogLikelihoods[i,state] <- sum(dbinom(Dt[i,],Dc[i,],p,log=TRUE))
      if(is.nan(rowDataLikelihoods[i,state])){
        warning("NaN!!!!!! rowDataLikelihoods[i,state]")
      }
      if(is.nan(rowDataLogLikelihoods[i,state])){
        warning("NaN!!!!!! rowDataLogLikelihoods[i,state]")
      }
      ## compute the product of Bernoulli distributions needed for 
      ## Pr(Z.i|mu)
      row.z.likelihoods[i,state] <- row.z.likelihood( z, Mu[i,] )
      if(doLog){
        qu[i,state] <- rowDataLogLikelihoods[i,state]+
          log(row.z.likelihoods[i,state])
      } else {
        qu[i,state] <- rowDataLikelihoods[i,state]*row.z.likelihoods[i,state]
      }
    }
  }
  result[["unNormalized"]] <- qu
  ##
  ##Normalizing:
  result[["normalized"]] <- normalizeQ(qu=qu, doLog=doLog,
                                       doEnforceNormal=doEnforceNormal)$qu
  normalizedUnLog <- result[["normalized"]]
  if(doLog)
    normalizedUnLog <- exp(normalizedUnLog)
  result[["normalizedUnLog"]] <- normalizedUnLog
  return(result)
}

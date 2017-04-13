Clomial.likelihood <-
function(Dc, Dt, Mu, P){
  ## Sums log-likelihoods assuming a binomial distribution for each Dt_[i,t]
  ## Where n = Dc_[i,t], k = Dt_[i,t], and p =  Mu_i*P[,t] 
  ## OUTPUT: ll is the log-likelihood, and llS is the matrix for all mutations
  ## and clones.
  result <- list()
  log.likelihood <- 0
  llS <- c()
  ##
  for(i in 1:nrow(Mu)){
    ##p <- Mu[i,]%*%P
    
    #p <- (Mu[i,]/2)%*%P ## /2 because of heterozygosity
    p <- (Mu[i,])%*%P
    
    p[p>1] <- 1
    ##
    lli <- dbinom(Dt[i,],Dc[i,],p,log=TRUE)
    log.likelihood <- log.likelihood + sum(lli)
    llS <- rbind(llS,lli)
  }
  rownames(llS) <- rownames(Mu)
  result[["ll"]] <- log.likelihood
  result[["llS"]] <- llS
  return(result)
}

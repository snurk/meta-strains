Clomial.generate.data <-
function(N, C, S, averageCoverage, mutFraction,doSample1Normal=FALSE, erroRate=0, doCheckDc=TRUE){
    ##
    ## doCheckDc: If TRUE, generating with be repeated until no row of Dc is all zeros
    ##^I.e. all loci should have positive coverage in at least one sample.
    ## QC:
  if(erroRate<0 | 1<erroRate )
    stop("erroRate is the probability of effective sequencing error,
 and should be in range [0,1]!")
  ## Randomly assign N*S*averageCoverage reads to each locus and sample
  Dc <- matrix(0,N,S)
  while(TRUE){
      rnums <- floor( runif( averageCoverage*N*S )*N*S )
      for( i in rnums ){
          Dc[i+1] <- Dc[i+1]+1
      }
      if(sum(rowSums(Dc)==0)==0 | !doCheckDc)
          break
  }
  ##
  ## Generate U such that the first column is 0s 
  ## and the other elements have a probability of mutFraction to be 1
  U <- matrix(runif(N*(C-1)),N,C-1) ## True genotypes of clones
  U[U<mutFraction] <- 1
  U[U != 1] <- 0
  U <- cbind( matrix(0,N,1), U ) ## the first column of U is 0s.
  ## Check to see that at least one clone for each locus contains 
	## a cancer allele. If not, re-assign values at that locus.
  for( i in 1:nrow(U) ){
    while( sum( U[i,] ) == 0 ){
      Ui <- matrix(runif(C-1),1,C-1)
      Ui[ Ui<mutFraction ] <- 1
      Ui[ Ui != 1 ] <- 0
      Ui <- cbind( 0, Ui )
      U[i,] <- Ui
    }
  }
  ##
  ## Generate P such that each column is a normalized vector of uniform draws
  PTrue <- matrix(runif(C*S),C,S)
  for( t in 1:S ){
    PTrue[,t] <- PTrue[,t]/sum(PTrue[,t]) ## normalizing the column.
  }
  if(doSample1Normal)
    PTrue[,1] <- c(1,rep(0,times=C-1))
  ##
  ## No sample can be too close to normal.
  P <- protect.P(P=PTrue)$protected
  ##
  ## Dt is the expected value of the binary distribution given U*P and Dc
  ##
  
  
  #Phi <- (U/2)%*%PTrue ## /2 because of heterozygosity 
  Phi <- (U)%*%PTrue 
  
  Dt <- floor(Dc * Phi) ## Only initiates Dt.
  ##
  ## Dt is random samples of the binary distribution given U*P and Dc
  if(TRUE ){
    for( i in 1:N ){
      for( t in 1:S ){
        noisyPhi <- Phi[i,t]*(1-2*erroRate) + erroRate
        ##^ phi*(1-e)+(1-phi)*e
        Dt[i,t] <- rbinom(1,Dc[i,t],noisyPhi)
      }
    }
  }
  if(doSample1Normal)
    Dt[,1] <- 0 ## The first sample is normal.
  ##
  data <- list()
  data[["Dc"]] <- Dc
  data[["Dt"]] <- Dt
  data[["U"]] <- U
  data[["PTrue"]] <- PTrue
  data[["Likelihood"]] <- Clomial.likelihood(Dc,Dt,U,PTrue)$ll
  data[["Phi"]] <- Phi
  ##
  return(data)
}

protect.P <-
function(P,maxNormal=0.9){
    ## Makes sure that each column of P (except the P[,1])
    ## has no more than maxNormal in the first row.
    result <- list()
    P2 <- P
    for(J in 2:ncol(P))
      if(P[1,J]>maxNormal){ ## This column needs protection.
        P2[,J] <- (1-maxNormal)/(nrow(P)-1)
        P2[1,J] <- maxNormal
      }##End if(P[1,J]>maxNormal).      
    result[["P"]] <- P
    result[["protected"]] <- P2
    return(result)
  }

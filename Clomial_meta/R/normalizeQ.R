normalizeQ <-
function(qu,doLog,doEnforceNormal=TRUE,scalar=50,doNormalize=TRUE){
    ## Normalizes q such that the sum of each row is 1.
    result <- list()
    N <- nrow(qu)
    C <- log(ncol(qu),base=2)
    newQ <- qu
    for( state in 1:(2^C) ){
      for( i in 1:N){
        s <- sum( qu[i,] )
        if( (!doLog & s==0) | (doLog & !is.finite(max(qu[i,]))) ){
          ## All probabilities are 0.
          if(!doEnforceNormal){
            qi <- rep(2^(-C),times=2^C)
          }else{
            ## The second half corresponds to normal and should be 0.
            qi <- rep(2^(-C),times=2^(C-1))*2
            qi <- c(qi,rep(0,times=2^(C-1)))
          }
          if(doLog){
            newQ[i,] <-log(qi)
          }else{
            newQ[i,] <-qi
          }
        }else
          if(doNormalize) { ## There is non-0 probability,
            if(!doLog){
              newQ[i,] <- qu[i,] / s
            }else{
              ## "Computation in log-space":
              ## Let x=c(log(A,B,C,D)),
              ## Ya=log(A/(A+B))=log(1/(A/A+B/A))=-log(A/A+B/A)
              ## B/A=exp(b)/exp(a)=exp(b-a)
              ## => Ya=-log(exp(a-a)+exp(b-a))=-log(sum(exp(x-a)))
              if(doLog){
                for(ind in 1:2^C){
                  if(is.finite(qu[i,ind])){
                    sum1 <- sum(exp((qu[i,]-qu[i,ind])+scalar)/exp(scalar))
                    newQ[i,ind] <- -log(sum1)
                  }else{
                    newQ[i,ind] <- -Inf
                  }
                }
              }
            }
          }
      }##End for(i).
    }##End for (state).
    result[["qu"]] <- newQ
    return(result)
  }

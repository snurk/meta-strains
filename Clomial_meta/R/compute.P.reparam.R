compute.P.reparam <-
function(C,S,N,Dt,Dc,ustates,qu,d2b,P,Mu=NULL,doEnforce=TRUE,trace=0,
           doShowRounds=FALSE,debug=FALSE,
           PTrue=NULL, doSilentOptim=TRUE,U=NULL,doDebug=FALSE,
           doTalk=TRUE){
    ## Computes the new P by first reparametrization and then BFGS-B;
    ## as explained in the paper.
    ## Mu and PTrue are included for debugging.  
    result <- list()
    if(is.null(rownames(P)))
      rownames(P) <- paste("C",1:nrow(P),sep="")
    if(is.null(colnames(P)))
      colnames(P) <- paste("S",1:ncol(P),sep="")
    optimumP <- P
    pre <- 0 ## total value of objective function before optimization
    post <- 0 ## after optimization
    bfgsRounds <- 0 ## max number of rounds objective function called by BFGS.
    for(J in 1:ncol(P)){ ## optimize each column (sample) independently.
      if(doDebug)
        message("J:",J)
      Pj <- P[,J]
      Wj <- P2W(Pj)
      checked <- check.Wj(Wj)
      ## The objective function and its gradient:
      ## Goal: maximize Phi.
      objective <- function(W)
        return(-Phi(Xj=Dt[,J],Rj=Dc[,J],Qu=qu,Pj=c(Pj[1],W2P(W))))
      gradient <- function(W)
        return(-dPhidW(Xj=Dt[,J],Rj=Dc[,J],Qu=qu,Wj=W)$dPhidWOut)
      ##
      if(doDebug|debug)
        browser()
      pre <- pre+objective(W=Wj)
      ##print(objective(W=Wj))
      ##
      update.Wj <- function(Wj,doDebug=FALSE){ ## helps in debugging.
        epsilon <- 10^(-7)
        resultPj <- list()
        optimized <- NULL
        ##
        lowers <- rep(epsilon,times=length(Wj))
        ## Optimizing:
        tried <- try(optimized <- optim(par=Wj,fn=objective,gr=gradient,
                                        lower=lowers,
                                        method="L-BFGS-B",
                                        control=list(trace=trace)),
                     silent=doSilentOptim)
        ##
        resultPj[["tried"]] <- tried
        resultPj[["optimized"]] <- optimized
        if(doDebug)
          browser()
        return(resultPj)
      }##End update.Wj <- function.
      ##
      updated <- update.Wj(Wj=Wj,doDebug=FALSE)
      ##^ the function helps in debugging.
      if(class(updated$tried)!="try-error") {
        ## otherwise do not update this column of P.
        rounds <- updated$optimized$counts["function"]
        bfgsRounds <- max(bfgsRounds,rounds)
        if(doShowRounds)
          print(paste("Optimized after rounds of: ",rounds))
        updatedWj <- updated$optimized$par
        optimumPj <- W2P(W=updatedWj)
        optimumPj <- c(1-sum(optimumPj),optimumPj)
        post <- post+objective(W=updatedWj)
        if(sum(optimumPj<0)!=0 & doEnforce) ##QC.
          stop("None appropriate update for P!")
        ##
        ## Updating:
        optimumP[,J] <- optimumPj
        plot.obj <- function(W,ind=1,st=0.01){
          ## plots the objective function by varying W
          ## st: the step
          result <- list()
          xs <- -100:100*st  
          ys <- c()
          one <- rep(0,length(W))
          one[ind] <- 1
          for(x1 in xs)
            ys[as.character(x1)] <- objective(W+one*x1)
          plot(x=xs,y=ys)
          result[["obj"]] <- ys ## the objective values around Wj
        }##End plot.obj <- function.
        if(debug)
          plot.obj(W=updatedWj)
      }
    }##End for(J in 1:ncol(P)).
    ##
    if(doTalk)
      print(paste(pre," -> ",post, ", ",round(100*(pre-post)/post,2),"%",
                  sep=""))
    optimumP <- protect.P(P=optimumP)$protected
    result[["P"]] <- P
    result[["post"]] <- post
    result[["pre"]] <- pre
    result[["optimumP"]] <- optimumP
    result[["bfgsRounds"]] <- bfgsRounds
    return(result)
  }

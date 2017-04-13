Clomial <-
function(Dt=NULL,Dc=NULL,DcDtFile=NULL,
                    C,doParal=FALSE,outPrefix=NULL,
                    binomTryNum=1000,maxIt=100,llCutoff=0.001,
                    jobNamePrefix="Bi",qstatWait=2,fitBinomJobFile=NULL,
                    jobShare=10,
                    ignoredSample=c(),
                    fliProb=0.05,conservative=TRUE, doTalk=FALSE){
  ## Runs many iterations to fit a binomial model to the data.
  ## It prepares the initializations and fits many models.
  ## -ignoredSample: The column index to be ignored.
  result <- list()
  binomTries <- list() ##The results of EM iterations to fit the binomials.
  ## qstat commands:
  qstatWaitCom <- paste("qstat |grep ",jobNamePrefix,"|grep qw|wc -l",sep="")
  qstatAll <- paste("qstat |grep ",jobNamePrefix,"|wc -l",sep="")
  ## Input check:
  if(is.null(Dt)|is.null(Dc)){ ## We expect to read it from file.
    if(is.null(DcDtFile))
      stop("Where should I read the counts from?")
    if(!file.exists(DcDtFile))
      stop(paste("No data file at:",DcDtFile))
    load(DcDtFile) ## Dt and Dc are loaded.
  }

  ## Ignoring some samples?
  if(length(ignoredSample)>0){
    if(doTalk)
      print(paste("These samples will be ignored:",ignoredSample))
    Dc <- Dc[,-ignoredSample]; Dt <- Dt[,-ignoredSample]
  }

  Dc <- as.matrix(Dc)
  Dt <- as.matrix(Dt)
  N <- nrow(Dc)
  S <- ncol(Dc)
  ## QC,
  if( (nrow(Dt)!=N) | (ncol(Dt)!=S) )
    stop("Dc and Dt should have the same dimension!")
  if( S<2 ){
    m1 <- "At least 2 samples are needed to run Clomial!"
    m2 <- "A normal sample can be simulated by adding"
    m3 <- "a column of 1s and 0s to Dc and Dt, accordingly."
   stop(paste(m1,m2,m3))
  }
  if(maxIt<3){
    stop("It is recommended to have at least 3 iterations (maxIt>2) !")
  }
  ##
  freq1 <- Dt/Dc
  if(0 %in% rowSums(Dc)){
      stop("Some rows of Dc are zero!")
  }
  for(J in 1:binomTryNum){
    if(doTalk)
      print(paste("--- Training the ",J,"th model...",sep=""))
    if(!doParal){
      ## Random initialization:
      random1 <- runif(n=N*(C-1),min=rowMins(freq1,na.rm=TRUE)*0.9,max=rowMaxs(freq1,na.rm=TRUE)*1.1)
      random1[random1>1] <- 1
      random1[random1<0] <- 0
      Mu <- matrix(random1,N,C-1)
      Mu <- cbind( matrix(0,N,1), Mu )
      rownames(Mu) <- rownames(Dc)
      colnames(Mu) <- paste("C",1:C,sep="")
      P <- matrix(runif(C*S),C,S)
      rownames(P) <- colnames(Mu)
      colnames(P) <- colnames(Dc)
      ## Normalize P:
      for( t in 1:S ){
        s <- sum(P[,t])
        P[,t] <- P[,t]/s
      }##End for.
      ##
      binomTries[[J]] <- try(Clomial.iterate(maxIt=maxIt, Dc=Dc, Dt=Dt,
                                             Mu=Mu, P=P,
                                             doSilentOptim=FALSE,
                                             llCutoff=llCutoff,
                                             computePFunction=compute.P.reparam,
                                             doTalk=doTalk,
                                             fliProb=fliProb,
                                             conservative=conservative))
    } else { ## Submit jobs to cluster
      dir.create(outPrefix,showWarnings=FALSE)
      tryName <- paste(jobNamePrefix,J,sep="")
      outPrefixJ <- paste(outPrefix,"/",tryName,sep="")
      while(TRUE){
        doSubmitNow <- TRUE
        if((J%%100)==0){
          totalJobsNum <- system(qstatAll,intern=TRUE)
          totalJobsNum <- as.numeric(totalJobsNum)
          waitingJobsNum <- system(qstatWaitCom,intern=TRUE)
          waitingJobsNum <- as.numeric(waitingJobsNum)
          if(waitingJobsNum > (1/2)*totalJobsNum+200){
            ## wait for some jobs to finish.
            doSubmitNow <- FALSE
            Sys.sleep(qstatWait)            
          }
        }
        if(doSubmitNow){
          print(J)
          args1 <- c(DcDtFile,0,C,llCutoff,outPrefixJ,fliProb,ignoredSample)
          s1 <- rsub(tempFolder=outPrefix,sourceFile=fitBinomJobFile,
                     namePrefix=tryName,wait=0.1,
                     arguments=args1,
                     memory="1.0G",jobShare=jobShare,
                     doSubmit=TRUE,doDeleteTempFile=FALSE,jobName=tryName,
                     ste=NULL,sto=NULL,
                     doTalk=FALSE,doWarn=FALSE)
          print(J)
          break 
        }
      }
    }## End else, a job was submitted.
  }##End for(J in 1:binomTryNum).
  if(doParal){
    ## These are the files where the outputs were saved in.
    binomTries <- paste(outPrefix,"/",jobNamePrefix,1:binomTryNum,".RData",
                        sep="")
  }
  if(!is.null(outPrefix)){
    save(binomTries,file=paste(outPrefix,".RData",sep=""))
  }
  result[["models"]] <- binomTries
  return(result)
}

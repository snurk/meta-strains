rsub <-
function(tempFolder="~/temp/", sourceFile,
                 namePrefix=basename(sourceFile),wait=60,
                 arguments=c(), memory="10.0G",jobShare=200,
                 sourceBashrcBy="~/.bashrc all",
                 doSubmit=TRUE,doDeleteTempFile=FALSE,jobName=NULL,
                 ste=NULL,sto=NULL,
                 doTalk=TRUE,doWarn=TRUE){
  ## Submitts a job to the cluster by qsub:
  ## -ste, sto: The files to save standard error, output in.
  startTime <- Sys.time()
  tempFolder <- paste(tempFolder,"/rsub/",basename(sourceFile),"/",sep="")
  dir.create(tempFolder,recursive=TRUE,showWarnings=FALSE)
  argumentsString <- paste(arguments,collapse="_")
  if(is.null(jobName))
    jobName <- paste(namePrefix,argumentsString,sep="_")
  tempFile <- paste(tempFolder,"/",jobName,"_",startTime,"_qsub.sh",sep="")
  tempFile <- gsub(x=tempFile,pattern=" ",replacement="_" )
  tempFile <- gsub(x=tempFile,pattern=":",replacement="-" )
  if(is.null(ste))
    ste <- paste(tempFile,".e",sep="") ## standard error
  if(is.null(sto))
    sto <- paste(tempFile,".o",sep="")
  ##
  args1 <- paste(arguments,collapse=" ")
  system(paste("echo 'source ",sourceBashrcBy," ' > ",tempFile,sep=""))
  system(paste("echo Rscript  ",sourceFile," ",args1 ," >> ",tempFile,sep=""))
  system(paste("chmod 777 ",tempFile))
  ## Making command:
  commanOptins <- paste("qsub -js ",jobShare,
                        " -l mem_free=",memory, " -l h_vmem=",
                        memory," -l mem_requested=", memory,
                        " -l longjob=TRUE -N ",
                        jobName," -e ",ste," -o ",sto,sep="")
  command <- paste(commanOptins," ",tempFile,sep="")
  if(doTalk)
    print(command)
  ##
  if(doSubmit)
    system(command)
  if(wait<1 & doWarn)
    warning("Jobs may overwrite if wait<1 ")
  Sys.sleep(wait) ## Wait between submitting jobs.
  if(doDeleteTempFile)
    system(paste("rm -rf ",tempFile,sep=""))
  return(list(tempFile=tempFile,jobName=jobName,command=command))
}

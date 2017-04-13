\name{Clomial.iterate}
\alias{Clomial.iterate}
\title{
	Runs EM iterations until convergence of the Clomial model.
}
\description{
  Given the data and the initial values for the model parameters,
  runs EM iterations until convergence of the Clomial model.
}
\usage{
Clomial.iterate(Dt, Dc, Mu, P, maxIt=100, U = NULL, PTrue = NULL,
llCutoff = 10^(-3), computePFunction = compute.P.reparam,
doSilentOptim = TRUE, doTalk = TRUE, doLog = TRUE, debug = FALSE,
noiseReductionRate = 0.01, fliProb=0.05,conservative=TRUE)
}
\arguments{
  \item{maxIt}{
	The maximum number of EM iterations.
}
  \item{Dt}{
	A matrix which contains the counts of the alternative allele
	where rows correspond to the genomic loci, and columns
	correspond to the samples.
}
  \item{Dc}{
	A matrix which contains the counts of the total number of mapped
	reads where rows correspond to the genomic loci, and columns
	correspond to the samples.
}
  \item{Mu}{
    The initial value for the Mu matrix which models the genotypes,
    where rows and columns correspond to genomic loci and clones, accordingly.
}
  \item{P}{
	The initial matrix of clonal frequency where rows and columns
	correspond to clones and samples, accordingly.
}
  \item{U}{
	The true value for \code{Mu}, used for debugging purposes only.
}
  \item{PTrue}{
	The true value for \code{P}, used for debugging purposes only.
}
  \item{llCutoff}{
	EM iterations stops if the relative improvement in the
	log-likelihood is not more than this threshold. 
}
  \item{computePFunction}{
    The function used for updating \code{P}.
    For advanced development use only.
}
  \item{doSilentOptim}{
	If TRUE, the optimization massages will not be reported.
}
  \item{doTalk}{
	If FALSE, the function will be run in silent mode.
}
  \item{doLog}{
    Highly recommended to set to TRUE. Then, the computations will
    be done in log space to avoid numerical issues. 
}
  \item{debug}{
	If TRUE, the debug mode will be turned on.
}
  \item{noiseReductionRate}{
	The noise will be reduce by this rate after each EM iteration.
}
	\item{fliProb}{
	A "flipping probability" used for noise injection which can be
	disabled when \code{fliProb=0}. After the first EM iteration,
	each entry of the matrix \code{Mu} such as \code{m} may change
	to \code{1-m} with this probability. This probability decreases
	on subsequent iterations. 
}
	\item{conservative}{
	Boolean where TRUE means noise will be injected only if
	likelihood is improved after an EM iteration, otherwise the
	original Mu matrix will be used for the next iteration. 
	For expert use only.
}
}
\details{
	Injecting noise can be done by assigning a positive value to
	\code{fliProb}, and can be disabled by \code{fliProb=0}.
	Noise injection is recommended for training models with a high
	number of clones (>4).
}

\value{
A list will be made with the following entries:
	\item{Qs }{The history of matrices containing the posterior
	  \code{Q} values.}
	\item{Ps }{The history of \code{P} matrices.}
	\item{Mus }{The history of \code{Mu} matrices.}
	\item{Mu }{The value of \code{Mu} after convergence.}
	\item{P }{The value of \code{P} after convergence.}
	\item{llCutoff }{The threshold used to decide convergence.}
	\item{LRatio }{The final relative improvement in the log
	  likelihood which lead to convergence.}
	\item{Likelihoods }{The history of log-likelihoods.}
	\item{fliProb }{The final value of \code{fliProb} used for noise
	  injection.}	
	\item{timeTaken }{An object of class ``\link{difftime}'' which
	  reports the total computational time for EM iterations.}
	\item{endTaken }{An object of class ``POSIXct''
	  (see \link{DateTimeClasses}) which reports the time EM
	  iterations finished.}
}
\references{
  Inferring clonal composition from multiple sections of a breast cancer,
  Zare et al., Submitted.
}
\author{
Habil Zare
}

\seealso{
\code{\link[Clomial:Clomial-package]{Clomial}}, 
\code{\link{Clomial}}, \code{\link{breastCancer}}
}

\examples{
set.seed(1)
## Getting data:
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt
freq1 <- Dt/Dc
N <- nrow(Dc)
S <- ncol(Dc)
Cnum <- 4 ## assumed number of clones.
## Random initialization:
random1 <- runif(n=N*(Cnum-1),min=rowMins(freq1)*0.9,max=rowMaxs(freq1)*1.1)
random1[random1>1] <- 1
random1[random1<0] <- 0
Mu <- matrix(random1,N,Cnum-1)
Mu <- cbind( matrix(0,N,1), Mu )
rownames(Mu) <- rownames(Dc)
colnames(Mu) <- paste("C",1:Cnum,sep="")
P <- matrix(runif(Cnum*S),Cnum,S)
rownames(P) <- colnames(Mu)
colnames(P) <- colnames(Dc)
## Normalizing P:
for( t in 1:S ){
	s <- sum(P[,t])
	P[,t] <- P[,t]/s
}##End for.

## Running EM: 
model1 <- Clomial.iterate(Dt=Dt, Dc=Dc, Mu=Mu, P=P)
print("Genotypes:")
round(model1$Mu)
print("Clone frequencies:")
model1$P
}
\keyword{ iteration }
\keyword{ documentation }
\keyword{ models }

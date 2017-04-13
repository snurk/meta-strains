\name{choose.best}
\alias{choose.best}
\title{
	Chooses the best trained Clomial model.
}
\description{
	Given the output of Clomial function, the likelihoods of all
	models are compared, and the best model is determined.
}
\usage{
choose.best(models, U = NULL, PTrue = NULL, compareTo = NULL, upto =
"All", doTalk=FALSE)
}
\arguments{
  \item{models}{
	The models trained by \code{Clomial} function.	
}
  \item{U}{
	The optional genotype matrix used for comparison. 
}
  \item{PTrue}{
	The optional clone frequency matrix used for comparison. 
}
  \item{compareTo}{
	The index of the model against which all other models are
	compared. Set to \code{NULL} to disable. 
}
  \item{upto}{
    The models with index less than this value are considered.
    Set to "All" to include every model.
}
	\item{doTalk}{
	If TRUE, information on number of analyzed models is reported.
}
}

\details{
  If \code{compareTo}, \code{U}, and \code{PTrue} are \code{NULL}
  no comparison will be done, and the function runs considerably faster. 
}
\value{
A list will be made with the following entries:
\item{err }{A list with 2 entries; err$P and err$U
  the vectors of clonal frequency errors, and genotype errors, accordingly.}
\item{Li }{A vector of the best obtained log-likelihood for each model.}
\item{bestInd }{The index of the best model in terms of log-likelihood.}
\item{comparison }{If \code{compareTo} is not \code{NULL}, the result
  of comparison with the corresponding model is reported.}
\item{bestModel }{The best model in terms of log-likelihood.}
\item{seconds }{A vector of the time taken,  in seconds, to train each model.}
}
\references{
  Inferring clonal composition from multiple sections of a breast cancer,
  Zare et al., Submitted.
}
\author{
Habil Zare
}
\note{
  When the number of assumed clones, \code{C}, is greater than 6,
  the comparison will be time taking because all possible permutations
  of clones should be considered. The running time will be slowed down
  by \code{C!}.
}

\seealso{
  \code{\link[Clomial:Clomial-package]{Clomial}},
  \code{\link{Clomial.likelihood}}, \code{\link{Clomial.iterate}}
}
\examples{
set.seed(4)
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt
ClomialResult <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=4,doParal=FALSE,binomTryNum=5)
chosen <- choose.best(models=ClomialResult$models)
M1 <- chosen$bestModel
print("Genotypes:")
round(M1$Mu)
print("Clone frequencies:")
M1$P
bestInd <- chosen$bestInd
plot(chosen$Li,ylab="Log-likelihood",type="l")
points(x=bestInd,y=chosen$Li[bestInd],col="red",pch=19)
}
\keyword{ documentation }
\keyword{ models }

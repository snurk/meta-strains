\name{Clomial1000}
\alias{Clomial1000}
\docType{data}
\title{
Pre-computed results of Clomial.}

\description{
Pre-computed results of Clomial function are provided for demo purposes.
It contains 1000 trained models on counts data from multiple samples of
a single primary breast cancer
obtained by deep, next-generation sequencing.}

\usage{data(Clomial1000)}
\format{\code{Clomial1000[["models"]]} is the list of trained models.}

\details{
Each model is the output of \code{Clomial.iterate()} function 
on the \code{breastCancer} data assuming there are 4 clones.}


\references{
  Inferring clonal composition from multiple sections of a breast cancer,
  Zare et al., Submitted.}

\seealso{
\code{\link[Clomial:Clomial-package]{Clomial}},
\code{\link{Clomial.iterate}}, \code{\link{choose.best}},
\code{\link{breastCancer}}
}

\examples{
data(Clomial1000)
chosen <- choose.best(models=Clomial1000$models)
M1 <- chosen$bestModel
print("Genotypes:")
round(M1$Mu)
print("Clone frequencies:")
M1$P
bestInd <- chosen$bestInd
plot(chosen$Li,ylab="Log-likelihood",type="l")
points(x=bestInd,y=chosen$Li[bestInd],col="red",pch=19)
}

\keyword{datasets}

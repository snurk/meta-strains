\name{compute.errors}
\alias{compute.errors}
\title{
	Computes the error of a Clomial model.
}
\description{
  Given the true genotype and frequency matrices,
  finds the permutation of genotypes matrix which best matches
  the true genotypes and returns the corresponding errors.}
      
\usage{
compute.errors(Mu, U, P, PTrue)
}
\arguments{
\item{Mu}{
  The matrix which models the genotypes, where rows and columns
  correspond to genomic loci and clones, accordingly.
}
\item{U}{
  The true genotype matrix defined similar to \code{Mu}.
}

\item{P}{
  The matrix of clonal frequency where rows and columns correspond
  to clones and samples, accordingly.
}
\item{PTrue}{
  The true clonal frequency matrix defined similar to \code{P}.
}
}

\details{
  Computing the error is useful for estimating the performance of
  inference on simulated, and for comparing different trained models.
  Genotype and frequency errors are defined as
  the normalized \code{l1}-error in
  reconstructing the genotype, and the clone frequency matrices,
  accordingly, where by normalized \code{l1}-error we mean the sum of absolute
  values of an error matrix divided by the size of the matrix.
}

\value{
  A list will be made with the following entries:
  \item{UError }{The \code{l1}-error of the genotype matrix
    normalized by the size of matrix.}
  \item{discretizedUError }{The \code{l1}-error of the rounded genotype
    matrix,
    i.e. the number of mismatching genotypes, normalized
    by the size of matrix}.
  \item{PErrorAbsolute }{The normalized \code{l1}-error of the clone frequency matrix.}
  \item{PErrorRelative }{Each entry of the error clone frequency
    matrix is normalized by the corresponding entry in
    \code{PTrue}, and then the normalized \code{l1} norm is computed.}
}
\references{
  Inferring clonal composition from multiple sections of a breast cancer,
  Zare et al., Submitted.
}
\author{
Habil Zare
}
\note{
  The use of \code{UError} and \code{PErrorAbsolute} is recommended.
  Computing the error is not feasible for more than 7 clones because
  the number of all possible permutations is factorial in the
  number of clones which grows super fast. Such input will trigger an
  error message.}

\seealso{
\code{\link[Clomial:Clomial-package]{Clomial}}
}
\examples{
set.seed(1)
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt
bics <- c()
ClomialResult <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=3,doParal=FALSE,binomTryNum=2)
model1 <- ClomialResult$models[[1]]
model2 <- ClomialResult$models[[2]]
## Comparing 2 trained models:
compute.errors(Mu=model1$Mu,U=model2$Mu,P=model1$P,PTrue=model2$P)
}

\keyword{ documentation }

\name{Clomial.likelihood}
\alias{Clomial.likelihood}
\title{
	Computes the complete data log-likelihood of a Clomial model.
}
\description{
	Computes the expected complete data log-likelihood of a Clomial
	model over all possible values of the hidden variables.
}
\usage{
Clomial.likelihood(Dc, Dt, Mu, P)
}
\arguments{
 \item{Dt}{
	A matrix which contains the counts of the alternative allele
	where rows correspond to the genomic loci, and columns
	correspond to the samples.
	
}
\item{Dc}{
  A matrix which contains the counts of the total number of mapped
  reads where rows correspond to the genomic loci, and columns correspond
  to the samples.
}
\item{Mu}{
  The matrix which models the genotypes, where rows and columns
  correspond to genomic loci and clones, accordingly.
}
\item{P}{
  The matrix of clonal frequency where rows and columns correspond to
  clones and samples, accordingly.
 }
}
\details{
	By assuming that the genomic loci and the samples are
	independent given the model parameters, the computation is
	simplified by first summing over the samples for a locus, and
	then summing over all the loci. This strategy avoids exploring
	the exponentially huge probability space. 
	
}
\value{
A list will be made with the following entries:
	\item{ll }{The expectation of complete log-likelihood over the
	  hidden variables.}
	\item{llS }{A vector of computed log-likelihoods at all loci.}
}
\references{
  Inferring clonal composition from multiple sections of a breast cancer,
  Zare et al., Submitted.
}
\author{
Habil Zare
}
\note{
The likelihood is computed assuming the heterozygosity is 2.
}

\seealso{
\code{\link[Clomial:Clomial-package]{Clomial}},
\code{\link{choose.best}},
\code{\link{compute.bic}}, \code{\link{breastCancer}}
}

\examples{
set.seed(1)
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt
ClomialResult <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=4,doParal=FALSE,binomTryNum=1)
model1 <- ClomialResult$models[[1]]
likelihood <- Clomial.likelihood(Dc=Dc, Dt=Dt, Mu=model1$Mu, P=model1$P)$ll
print(likelihood)
}

\keyword{ documentation }
\keyword{ models }

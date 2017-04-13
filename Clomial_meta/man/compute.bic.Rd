\name{compute.bic}
\alias{compute.bic}
\title{
	Computes BIC for a Clomial model.
}
\description{
	Computes the Bayesian Information Criterion (BIC) for a Clomial
	model, which might be useful to estimate the number of clones. A
	"significantly" smaller BIC is usually interpreted as a better
	fit to the data.
}
\usage{
	compute.bic(Dc, Dt, Mu, P)
}
\arguments{
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
	The matrix which models the genotypes, where rows and columns
	correspond to genomic loci and clones, accordingly.
}
  \item{P}{
	The matrix of clonal frequency where rows and columns correspond
	to clones and samples, accordingly.
}
}

\details{
	The Bayesian Information Criterion (BIC) for a model is computed
	by subtracting the expected log-likelihood times 2, from the
	number of free parameters of the model times logarithm of the
	total number of observations. For a Clomial model, we have
	\code{BIC = (NC+SC-S)log(sum(Dc))-2L}, where \code{L} is the
	likelihood, \code{N} is the number of genomic loci, \code{C} is
	the assumed number of clones, \code{S} is the number of samples,
	and \code{sum(Dc)} is the total number of observed reads.   
}
\value{
A list will be made with the following entries:
	\item{bic }{The BIC value.}
	\item{aic }{The AIC value.}
	\item{obsNum }{The total number of observed reades.}
}
\references{
  Inferring clonal composition from multiple sections of a breast cancer,
  Zare et al., Submitted.
}
\author{
Habil Zare
}
\note{
Theoretically, a method such as the Bayesian information criterion (BIC)
	or the Akaike information criterion (AIC) may be applied to
	estimate the number of clones. However, in practice, the outcome
	of such approaches should be interpreted with great caution
	because some of the underlying assumptions of the statistical
	analysis may not be necessarily true for a given model. For
	example, while a "small" improvement in the BIC is generally
	considered as a sign to stop making the model more complicated,
	making such decisions is very objective, and requires relying on
	thresholds with little statistical basis.
}

\seealso{
\code{\link[Clomial:Clomial-package]{Clomial}}
}
\examples{
set.seed(1)
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt
bics <- c()
Clomial3 <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=3,doParal=FALSE,binomTryNum=1)
model3 <- Clomial3$models[[1]]
bics[3] <- compute.bic(Dc=Dc,Dt=Dt, Mu=model3$Mu, P=model3$P)$bic
Clomial4 <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=4,doParal=FALSE,binomTryNum=1)
model4 <- Clomial4$models[[1]]
bics[4] <- compute.bic(Dc=Dc,Dt=Dt, Mu=model4$Mu, P=model4$P)$bic
print(bics) ## 4 is a better estimate for the number of clones.
}

\keyword{ documentation }

\name{Clomial.generate.data}
\alias{Clomial.generate.data}
\title{
  Generates simulated data to test performance of Clomial algorithm.
}
\description{
  Data sets are simulated based on binomial distribution using
  random parameters for the model. The accuracy of the EM procedure
  can be estimated by comparing the inferred parameters vs. the known
  ones which were used to generate the data.
}
\usage{
Clomial.generate.data(N, C, S, averageCoverage, mutFraction,
doSample1Normal = FALSE,erroRate=0,doCheckDc=TRUE)
}
\arguments{
  \item{N}{
	The number of genomic loci.
}
  \item{C}{
	The number of clones.
}
  \item{S}{
	The number of samples.
}
  \item{averageCoverage}{
	The average coverage over each loci, each sample.
}
  \item{mutFraction}{
    Should be in range 0-1. Each loci in every sample can be mutated
    with this probability.
}
  \item{doSample1Normal}{
    If TRUE, no contamination with the tumor content is allowed for
    the normal sample. I.e. the first column of the generated \code{P}
    matrix will start with 1, and the rest of its entries will be equal to 0. 
}
  \item{erroRate}{
    The sequencing noise can be simulated by assigning a positive value
    to this parameter, which is the probability of reading a normal
    allele as the alternative allele, and vica versa.
}
  \item{doCheckDc}{
  If TRUE, generating with be repeated until no row of Dc is all zeros to guarantee
  all loci have positive coverage in at least one sample.
}

}
\details{
See the reference below for details.
}

\value{
A list will be made with the following entries:
	\item{Dc }{A matrix of simulated coverage for all loci and samples.}
	\item{Dt }{A matrix of alternative allele counts for all loci
	  and samples.}
	\item{Ptrue }{The true clone frequency matrix used for
	  generating the data.}
	\item{U }{The true genotype matrix used for generating the data.}
	\item{Likelihood }{The log-likelihood of the model with the
	  true parameters.}
	\item{Phi}{The matrix of the second parameters of the binomial
	  distributions; each entry is the probability that a read
	  contains the variant allele at a locus in a sample.}
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
  \code{\link{Clomial.likelihood}}
}
\examples{
set.seed(1)
simulated <- Clomial.generate.data(N=20, C=4, S=10,
  averageCoverage=1000, mutFraction=0.1)
simulated$Dc
}
\keyword{ documentation }
\keyword{ models }

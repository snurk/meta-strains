\name{Clomial}
\alias{Clomial}
\title{
Fits several binomial models to data from multiple samples of a single tumor.
}
\description{
Using EM, trains several models using different initial values to escape
from local optima. The best one in terms of the likelihood can be later
chosen by choose.best() function.
}
\usage{
Clomial(Dt = NULL, Dc = NULL, DcDtFile = NULL, C, doParal=FALSE,
outPrefix = NULL, binomTryNum = 1000, maxIt = 100, llCutoff = 0.001,
jobNamePrefix = "Bi", qstatWait = 2, fitBinomJobFile = NULL,
jobShare = 10, ignoredSample = c(),  fliProb=0.05, conservative=TRUE,
doTalk=FALSE)

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
  \item{DcDtFile}{
	A file from which the data can optionally be loaded. It should
	contain the matrices Dc and Dt.
}
  \item{C}{
	The assumed number of clones.
}
  \item{doParal}{
	Boolean where TRUE means, in Linux, models with different
	initialization are trained in parallel on a cluster using qsub.
}
  \item{outPrefix}{
	A prefix for the path to save the results.
}
  \item{binomTryNum}{
	The number of models trained using different initialization.
}
  \item{maxIt}{
	The maximum number of EM iterations.
}
  \item{llCutoff}{
	EM iterations stops if the relative improvement in the
	log-likelihood is not more than this threshold. 
}
  \item{jobNamePrefix}{
	If run in parallel, this prefix will be used to name the jobs on
	the cluster.
}
  \item{qstatWait}{
	The waiting time between qstat commands to assess the number of
	running and waiting jobs.
}
  \item{fitBinomJobFile}{
	If run in parallel, this is the script which loads data, trains
	a model using a random initialization, and saves the results.
}
  \item{jobShare}{
	If run in parallel, the job_share option of qsub determines the
	priority of jobs over other submitted jobs.
}
  \item{ignoredSample}{
	A vector of indices of samples which will be ignored in
	training. Used by experts only to measure the stability of the
	results.
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
	original Mu matrix will be used for the next iteration. For
	expert use only.
}
	\item{doTalk}{
	If TRUE, information on the EM optimization iterations is reported.
}
}

\details{
	The likelihood of the model, given the hidden variables and the
	parameters, can be computed based on a combination of binomial
	distributions. In each EM iteration, the likelihood is
	increased, however, due to presence of local optima, several
	models should be tried using different random
	initialization. For higher number of assumed clones, \code{C},
	the parameter binomTryNum should be increased because the
	dimension of the search space grows linearly with \code{C}.
}

\value{
Returns a list containing the entry called \code{models},
which is a list of the length equal to binomTryNum where each element is
a trained model.
For each trained model, \code{Mu} models the matrix of genotypes, where
rows and columns correspond to genomic loci and clones,
accordingly. Also, \code{P} is the matrix of clonal frequency where rows
and columns correspond to clones and samples, accordingly.
The first column of \code{P} corresponds to the normal clone. 
The history of \code{Mu}, \code{P}, and the log-likelihood over
iterations is saved in lists \code{Ps}, \code{Mus}, and
\code{Likelihoods}, accordingly.
}
\references{
  Inferring clonal composition from multiple sections of a breast cancer,
  Zare et al., Submitted.
}
\author{
Habil Zare
}
\note{
	The parallel mode works only in Linux, and when qsub and qstat
	commands are available on a cluster.
}

\seealso{
\code{\link[Clomial:Clomial-package]{Clomial}},
\code{\link{choose.best}}, \code{\link{Clomial.iterate}},
\code{\link{compute.bic}}, \code{\link{breastCancer}}
}
\examples{
set.seed(1)
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt
ClomialResult <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=4,binomTryNum=2)
chosen <- choose.best(models=ClomialResult$models)
M1 <- chosen$bestModel
print("Genotypes:")
round(M1$Mu)
print("Clone frequencies:")
M1$P
}

\keyword{ iteration }
\keyword{ documentation }
\keyword{ models }

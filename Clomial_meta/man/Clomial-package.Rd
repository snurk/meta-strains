\name{Clomial-package}
\alias{Clomial-package}
\docType{package}
\title{
Fits a binomial model to data from multiple samples of a single tumor to
infer its clonal decomposition.
}
\description{
Clomial fits binomial distributions to counts obtained from Next Gen
Sequencing data of multiple samples of the same tumor. The trained
parameters can be interpreted to infer the clonal structure of the
tumor.
}
\details{
\tabular{ll}{
Package: \tab Clomial\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-02-11\cr
License: \tab GPL (>= 2)\cr
}
The main function is Clomial() which requires 2 matrices Dt and Dc among
its inputs. They contain the counts of the alternative allele, and the
total number of processed reads, accordingly. Their rows correspond to
the genomic loci, and their columns correspond to the samples. Several
models should be trained using different initial values to escape from
local optima, and the best one in terms of the likelihood can be chosen
by choose.best() function.
}
\author{
Habil Zare and Alex Hu

Maintainer: Habil Zare <zare@u.washington.edu>
}
\references{
  Inferring clonal composition from multiple sections of a breast cancer,
  Zare et al., PLoS Computational Biology 10.7 (2014): e1003703.
}
\keyword{ package }
\keyword{ datasets }
\keyword{ iteration }
\keyword{ documentation }
\keyword{ models }

\seealso{
\code{\link{Clomial}}, \code{\link{choose.best}},
\code{\link{Clomial.iterate}},  \code{\link{Clomial.likelihood}},
\code{\link{compute.bic}}, \code{\link{breastCancer}}
}
\examples{
set.seed(1)
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt
ClomialResult <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=4,doParal=FALSE,binomTryNum=2)
chosen <- choose.best(models=ClomialResult$models)
M1 <- chosen$bestModel
print("Genotypes:")
print(round(M1$Mu))
print("Clone frequencies:")
print(M1$P)
}

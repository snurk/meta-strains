\name{breastCancer}
\alias{breastCancer}
\docType{data}
\title{
Breast cancer data for clonal decomposition.}

\description{
Counts data from multiple samples of a single primary breast cancer obtained by deep, next-generation sequencing. 
The file is consist of two matrices \code{Dt} and \code{Dc} which contain the counts of the alternative alleles, and the total number of counts on each genomic loci for every tumor samples, accordingly.}


\usage{data(breastCancer)}
\format{A list containing 2 matrices.}

\details{
Each matrix contains counts of reads mapped to 17 genomic loci for 12 tumor samples where
the column A5-2 corresponds to the normal sample.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\references{
Inferring clonal composition from multiple sections of a breast cancer, Zare et al., Submitted.}

\seealso{
\code{\link[Clomial:Clomial-package]{Clomial}}
}
\examples{
data(breastCancer)
breastCancer$Dt
}
\keyword{datasets}

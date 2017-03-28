library('matrixStats', lib.loc="~/R/libs")
library('Clomial', lib.loc="~/R/libs")

args <- commandArgs(trailingOnly=TRUE)

R <- read.table(args[1])
X <- read.table(args[2])

ClomialResult <- Clomial(
  Dc = R, Dt = X,
  maxIt = 500,
  C = 4,
  binomTryNum = 10
)

chosen <- choose.best(models = ClomialResult$models)
M1 <- chosen$bestModel

write.table(round(M1$P, 2), args[3])

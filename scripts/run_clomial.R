library('matrixStats', lib.loc = "~/R/libs")
library('Clomial', lib.loc = "~/R/libs")

args <- commandArgs(trailingOnly = TRUE)

R <- read.table(args[1])
X <- read.table(args[2])
cur_C <- as.numeric(args[3])

meta.log.likelihood <- function(Dc, Dt, Mu, P) {
  log.likelihood <- 0
  for (i in 1:nrow(Mu)) {
    p <- (Mu[i,]) %*% P
    p[p > 1] <- 1
    lli <-
      dbinom(as.numeric(Dt[i,]), as.numeric(Dc[i,]), p, log = TRUE)
    log.likelihood <- log.likelihood + sum(lli)
  }
  return(log.likelihood)
}


meta.bic_and_ll <- function (Dc, Dt, Mu, P) {
  result <- list()
  ll <- meta.log.likelihood(Dc = Dc,
                            Dt = Dt,
                            Mu = Mu,
                            P = P)
  result[["ll"]] <- ll
  N <- nrow(Mu)
  C <- ncol(Mu)
  M <- ncol(P)
  obsNum <- sum(Dc)
  bic <- -2 * ll + (N * C + M * (C - 1)) * log(obsNum)
  result[["bic"]] <- bic
  return(result)
}

clomial_result <- Clomial(Dc = R,
                          Dt = X,
                          C = cur_C,
                          doTalk = TRUE)
chosen <- choose.best(models = clomial_result$models, doTalk = TRUE)
best_model <- chosen$bestModel
res <- meta.bic_and_ll(
  Dc = R,
  Dt = X,
  Mu = best_model$Mu,
  P = best_model$P
)
write.table(round(best_model$P, 2), paste("clomial_results/frequencies_", cur_C, ".txt", sep=""))
write.table(best_model$Mu, paste("clomial_results/genotypes_", cur_C, ".txt", sep=""))
write(res$bic, paste("clomial_results/BIC_", cur_C, ".txt", sep=""))

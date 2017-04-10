library('matrixStats', lib.loc="~/R/libs")
library('Clomial', lib.loc="~/R/libs")

args <- commandArgs(trailingOnly=TRUE)

R <- read.table(args[1])
X <- read.table(args[2])


meta.log.likelihood <- function(Dc, Dt, Mu, P) {
  log.likelihood <- 0
  for (i in 1:nrow(Mu)) {
    p <- (Mu[i, ]) %*% P
    p[p > 1] <- 1
    lli <- dbinom(as.numeric(Dt[i, ]), as.numeric(Dc[i, ]), p, log = TRUE)
    log.likelihood <- log.likelihood + sum(lli)
  }
  return(log.likelihood)
}


meta.bic_and_ll <- function (Dc, Dt, Mu, P) {
  result <- list()
  ll <- meta.log.likelihood(Dc=Dc, Dt=Dt, Mu=Mu, P=P)
  result[["ll"]] <- ll
  N <- nrow(Mu)
  C <- ncol(Mu)
  M <- ncol(P)
  obsNum <- sum(Dc)
  bic <- -2 * ll + (N * C + M * (C - 1)) * log(obsNum)
  result[["bic"]] <- bic
  return(result)
}


plot_bic_and_ll <- function() {
  bics <- c()
  log_likelihoods <- c()
  print("calculating likelihood and BIC for different number of strains...")
  for (i in c(3:10)) {
    print(i)
    
    clomial_result <- Clomial(Dc=R, Dt=X, C=i, doTalk=TRUE)
    chosen <- choose.best(models = clomial_result$models, doTalk=TRUE)
    best_model <- chosen$bestModel
    res <- meta.bic_and_ll(Dc=R, Dt=X, Mu=best_model$Mu, P=best_model$P)
    bics[i] <- res$bic
    log_likelihoods[i] <- res$ll
    
    png(filename="BIC.png")
    plot(bics, xlab="number of clones", ylab="BIC",
         col="red", pch=20, cex=2)
    dev.off()
    
    png(filename="LL.png")
    plot(log_likelihoods, xlab="number of clones", ylab="log(likelihood)",
         col="blue", pch=20, cex=2)
    dev.off()
  }
}

plot_bic_and_ll()

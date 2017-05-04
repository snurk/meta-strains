args <- commandArgs(trailingOnly = TRUE)

min_C <- as.numeric(args[1])
max_C <- as.numeric(args[2])

BICs <- c()
min_BIC <- Inf
res_C <- -1

for (cur_C in c(min_C:max_C)) {
  file_name <- paste("clomial_results/BIC_", cur_C, ".txt", sep = "")
  cur_BIC <- scan(file = file_name, what = double(), n = 1)
  BICs <- c(BICs, cur_BIC)
  print(cur_C)
  print(cur_BIC)
  if (cur_BIC < min_BIC) {
    min_BIC <- cur_BIC
    res_C <- cur_C
  }
}

print(paste("min BIC is in ", res_C, sep = ""))

png(filename=args[3])
plot(c(min_C:max_C), BICs,
     xlab = "number of clones",
     ylab = "BIC",
     pch = 20, col = "red", cex = 3,
     main = paste("min BIC is in ", res_C, sep = ""))
dev.off()
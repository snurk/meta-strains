setwd("/home/margarita/AU/meta_strains_project/benchmark/")
library(reshape2)
library(ggplot2)


JSD <- function(p, q){
  m <- 0.5 * (p + q)
  0.5 * (sum(p[p > 0] * log(p[p > 0]/ m[p > 0])) + sum(q[q > 0]* log(q[q > 0]/ m[q > 0])))
}

JSDmatr <- function(p, q){
  res <- rep(NA, ncol(p))
  for (i in 1:length(res)) {
    if (sum(p[, i] > 0) > 0)
    {
      res[i] <- JSD(p[, i], q[, i])
    }
  }
  res
}


x4 <- read.csv("4strains.prof", header = TRUE, sep = " ")
melted <- melt(x4)

#ggplot(melted, aes(x = method, y = value, fill = strain)) + 
#  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ variable) + labs(y="frequencies")

ideal4 <- x4[1:4, 2:12]
we4 <- x4[5:8, 2:12]
constr4 <- x4[9:12, 2:12]

print(ideal4)
print(we4)
print(constr4)


constrains4 <- print(JSDmatr(constr4, ideal4))
clomial4 <- print(JSDmatr(we4, ideal4))


test4 <- data.frame(samples = 1:11, constrains = constrains4, we = clomial4)

test24 <- melt(test4, id.var='samples')
test24 <- na.omit(test24)

names(test24) <- c("samples", "method", "value")
print(test24)

p4 <- ggplot(test24, aes(x=samples, y=value, color=method)) + 
  geom_line(size = 1.2) + geom_point(size = 3) + ggtitle("4 strains") +
  theme(plot.title = element_text(hjust = 0.5)) +
  
  xlab("sample") + ylab("Jensen–Shannon divergence") +
  scale_x_continuous(limits=c(1, 11), breaks = seq(1, 11, 1)) + 
  scale_y_continuous(limits=c(0, 0.4))

plot(p4)
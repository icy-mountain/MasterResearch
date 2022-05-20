rm(list = ls(all.names = TRUE))
install.packages('Rsolnp')
library(Rsolnp)
source("./utilities.R")
tab1_freq <- c(2, 3, 2, 1,
               1, 8, 6, 5,
               2, 4, 17, 5,
               1, 3, 3, 10)
tab3_freq <- c(248, 36, 5, 10,
               36, 49, 23, 15,
               4, 11, 13, 9,
               1, 1, 1, 9)
## "..." is solution of nsolnp bug
ENiMSConstrFunc <- function(p, ...) {
  rows <- CountRow(p)
  ENiMS_0Sum_Constr <- c(sum(p) - 1)
  #delta <- sum(ExtractW1Area(p, 1)) / (sum(ExtractW2Area(p, 1)))
  delta <- SumAllW1Area(p) / SumAllW2Area(p)
  for (i in 1:rows) {
    W1 <- sum(ExtractW1Area(p, i))
    W2 <- sum(ExtractW2Area(p, i))
    ENiMS_0Sum_Constr <- append(ENiMS_0Sum_Constr, W1 - delta * W2 )
  }
  return (ENiMS_0Sum_Constr)
}
ENiMSModel <- function(freq) {
  p0 <- freq / sum(freq)
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  eqB <- rep(0, rows + 1)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = ENiMSConstrFunc, eqB = eqB,
                      LB = lowerBound, freq = freq)
  solnpResult$df <- rows - 1
  return(solnpResult)
}
DisplayResult <- function(freq) {
  result <- ENiMSModel(freq = freq)
  phat <- result$pars
  delta <- sum(ExtractW1Area(phat, 1)) / sum(ExtractW2Area(phat, 1))
  result$modelParams <- delta
  mhat <- result$pars * sum(freq)
  result$G2 <- CalcG2(freq, mhat)
  result$pValue <- pchisq(q=result$G2, df=result$df, lower.tail=FALSE)
  result$AICp <- CalcAICplus(freq, result)
  print(sprintf("df:%s", result$df))
  print(sprintf("G2:%s", result$G2))
  print(sprintf("pValue:%s", result$pValue))
  print(sprintf("AICp:%s", result$AICp))
  return(result)
}
system.time(ENiMS_tab1_result <- DisplayResult(tab1_freq))
system.time(ENiMS_tab3_result <- DisplayResult(tab3_freq))
ENiMS_tab1_result$modelParams
ENiMS_tab3_result$modelParams
#######
rows <- CountRow(tab1_freq)
p0 <- rep(1/length(tab1_freq), rows*rows)
p0 <- tab1_freq / sum(tab1_freq)
tab1_freq[1:4] / sum(tab1_freq[1:4])
sumW1 <- 0
sumW2 <- 0
for (i in 1:rows) {
  sumW1 <- sumW1 + sum(ExtractW1Area(p0, i))
  sumW2 <- sumW2 + sum(ExtractW2Area(p0, i))
}
sumW1
sumW2

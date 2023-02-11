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
NiMSConstrFunc <- function(p, ...) {
  rows <- CountRow(p)
  NiMS_0Sum_Constr <- c(sum(p) - 1)
  for (i in 1:rows) {
    W1 <- sum(ExtractW1Area(p, i))
    W2 <- sum(ExtractW2Area(p, i))
    NiMS_0Sum_Constr <- append(NiMS_0Sum_Constr, W1 - W2)
  }
  return (NiMS_0Sum_Constr)
}
NiMSModel <- function(freq) {
  p0 <- rep(1/length(freq), length(freq))
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = NiMSConstrFunc,
                     eqB = rep(0, rows+1),
                     LB = lowerBound, freq = freq)
  solnpResult$df <- rows
  return(solnpResult)
}
DisplayNiMSResult <- function(freq) {
  result <- NiMSModel(freq = freq)
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
NiMS_tab1_result <- DisplayNiMSResult(tab1_freq)
NiMS_tab3_result <- DisplayNiMSResult(tab3_freq)

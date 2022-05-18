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
freq <- tab1_freq
### "..." is solution of nsolnp bug
NiMSConstrFunc <- function(p, ...) {
  rows <- CountRow(p)
  NiMS_0Sum_Constr <- c(sum(p) - 1)
  for (i in 1:rows) {
    W1 <- sum(CalcW1Area(p, i))
    W2 <- sum(CalcW2Area(p, i))
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
DisplayResult <- function() {
  result <- NiMSModel(freq = freq)
  df <- result$df
  mhat <- result$pars * sum(freq)
  G2 <- CalcG2(freq, mhat)
  pValue <- pchisq(q=G2, df=result$df, lower.tail=FALSE)
  AICp <- CalcAICplus(freq, result)
  print(sprintf("df:%s", df))
  print(sprintf("G2:%s", G2))
  print(sprintf("pValue:%s", pValue))
  print(sprintf("AICp:%s", AICp))
  return(result)
}
result <- DisplayResult()

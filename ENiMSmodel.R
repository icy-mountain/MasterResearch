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
ENiMSConstrFunc <- function(p, ...) {
  rows <- CountRow(p)
  ENiMS_0Sum_Constr <- c(sum(p) - 1)
  delta <- sum(ExtractW1Area(p, 1)) / (sum(ExtractW2Area(p, 1)))
  for (i in 2:(rows)) {
    W1 <- sum(ExtractW1Area(p, i))
    W2 <- sum(ExtractW2Area(p, i))
    ENiMS_0Sum_Constr <- append(ENiMS_0Sum_Constr,  W1 - delta * W2)
  }
  return (ENiMS_0Sum_Constr)
}
ENiMSModel <- function(freq) {
  p0 <- rep(1/length(freq), length(freq))
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  eqB <- rep(0, rows)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = ENiMSConstrFunc, eqB = eqB,
                      LB = lowerBound, freq = freq)
  solnpResult$df <- rows - 1
  return(solnpResult)
}
DisplayGNiMSResult <- function(freq) {
  result <- ENiMSModel(freq = freq)
  pHat <- result$pars
  delta <- sum(ExtractW1Area(pHat, 1)) / sum(ExtractW2Area(pHat, 1))
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
system.time(ENiMS_tab1_result <- DisplayGNiMSResult(tab1_freq))
system.time(ENiMS_tab3_result <- DisplayGNiMSResult(tab3_freq))
ENiMS_tab1_result$modelParams
ENiMS_tab3_result$modelParams

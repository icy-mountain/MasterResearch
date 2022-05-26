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
GNiMSConstrFunc <- function(p, ...) {
  rows <- CountRow(p)
  GNiMS_0Sum_Constr <- c(sum(p) - 1)
  delta <- sum(ExtractW1Area(p, 1)) / (sum(ExtractW2Area(p, 1)))
  numerator <- sum(ExtractW1Area(p, 2)) / sum(ExtractW2Area(p, 2))
  phi <- numerator / delta
  for (i in 3:rows) {
    W1 <- sum(ExtractW1Area(p, i))
    W2 <- sum(ExtractW2Area(p, i))
    GNiMS_0Sum_Constr <- append(GNiMS_0Sum_Constr, 
                                W1 - delta * phi^(i-1) * W2)
  }
  return (GNiMS_0Sum_Constr)
}
GNiMSModel <- function(freq) {
  p0 <- rep(1/length(freq), length(freq))
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  eqB <- rep(0, rows - 1)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = GNiMSConstrFunc, eqB = eqB,
                       LB = lowerBound, freq = freq)
  solnpResult$df <- rows - 2
  return(solnpResult)
}
DisplayGNiMSResult <- function(freq) {
  result <- GNiMSModel(freq)
  pHat <- result$pars
  delta <- sum(ExtractW1Area(pHat, 1)) / sum(ExtractW2Area(pHat, 1))
  numerator <- sum(ExtractW1Area(pHat, 2)) / sum(ExtractW2Area(pHat, 2)) 
  phi <- numerator / delta
  result$modelParams <- c(delta, phi)
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
system.time(GNiMS_tab1_result <- DisplayGNiMSResult(tab1_freq))
system.time(GNiMS_tab3_result <- DisplayGNiMSResult(tab3_freq))
GNiMS_tab1_result$modelParams
GNiMS_tab3_result$modelParams

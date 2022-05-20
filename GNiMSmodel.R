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
  delta <- sum(ExtractW1Area(p, 1)) / sum(ExtractW2Area(p, 1))
  numerator <- sum(ExtractW1Area(p, 2)) / sum(ExtractW2Area(p, 2))
  phai <- numerator / delta
  for (i in 1:rows) {
    W1 <- sum(ExtractW1Area(p, i))
    W2 <- sum(ExtractW2Area(p, i))
    GNiMS_0Sum_Constr <- append(GNiMS_0Sum_Constr, 
                                W1 - delta * phai^(i-1) * W2)
  }
  return (GNiMS_0Sum_Constr)
}
GNiMSModel <- function(freq) {
  p0 <- freq / sum(freq)
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  eqB <- rep(0, rows + 1)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = GNiMSConstrFunc, eqB = eqB,
                       LB = lowerBound, freq = freq)
  solnpResult$df <- rows - 2
  return(solnpResult)
}
DisplayResult <- function(freq) {
  result <- GNiMSModel(freq)
  phat <- result$pars
  delta <- sum(ExtractW1Area(phat, 1)) / sum(ExtractW2Area(phat, 1))
  numerator <- sum(ExtractW1Area(phat, 2)) / sum(ExtractW2Area(phat, 2)) 
  phai <- numerator / delta
  result$modelParams <- c(delta, phai)
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
system.time(GNiMS_tab1_result <- DisplayResult(tab1_freq))
system.time(GNiMS_tab3_result <- DisplayResult(tab3_freq))
GNiMS_tab1_result$modelParams
GNiMS_tab3_result$modelParams

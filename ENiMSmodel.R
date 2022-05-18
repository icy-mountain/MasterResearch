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
## "..." is solution of nsolnp bug
ENiMSConstrFunc <- function(p, ...) {
  rows <- CountRow(p)
  ENiMS_0Sum_Constr <- c(sum(p) - 1)
  for (i in 1:rows) {
    W1 <- sum(CalcW1Area(p, i))
    W2 <- sum(CalcW2Area(p, i))
    ENiMS_0Sum_Constr <- append(ENiMS_0Sum_Constr, W1 / (W2 + 1e-15))
  }
  return (ENiMS_0Sum_Constr)
}
ENiMSModel <- function(delta, freq) {
  p0 <- rep(1/length(freq), length(freq))
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  eqB <- c(0)
  eqB <- append(eqB, rep(delta, rows))
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = ENiMSConstrFunc, eqB = eqB,
                      LB = lowerBound, freq = freq)
  solnpResult$df <- rows - 1
  return(solnpResult)
}
forDeltaOptim <- function(delta, freq, output=FALSE) {
  if (output == TRUE)
    solnpResult <- ENiMSModel(delta, freq)
  if (output == FALSE) {
    sink(nullfile())
    solnpResult <- ENiMSModel(delta, freq)
    sink()
  }
  optimizedFuncValue <- solnpResult$value[length(solnpResult$value)]
  return(optimizedFuncValue)
}
DisplayResult <- function() {
  optimValues <- optimize(forDeltaOptim, interval = c(0, 10), freq = freq)
  optimDelta <- optimValues$minimum
  result <- ENiMSModel(delta = optimDelta,freq = freq)
  result$modelParams <- optimDelta
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
system.time(result <- DisplayResult())
result$modelParams

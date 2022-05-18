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
GNiMSConstrFunc <- function(p, phai, ...) {
  rows <- CountRow(p)
  GNiMS_0Sum_Constr <- c(sum(p) - 1)
  for (i in 1:rows) {
    W1 <- sum(CalcW1Area(p, i))
    W2 <- sum(CalcW2Area(p, i))
    GNiMS_0Sum_Constr <- append(GNiMS_0Sum_Constr, 
                                W1 / (phai^(i-1) * W2 + 1e-15))
  }
  return (GNiMS_0Sum_Constr)
}
GNiMSModel <- function(delta, phai, freq) {
  p0 <- rep(1/length(freq), length(freq))
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  eqB <- c(0)
  eqB <- append(eqB, rep(delta, rows))
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = GNiMSConstrFunc, eqB = eqB,
                       LB = lowerBound, freq = freq, phai = phai)
  solnpResult$df <- rows - 2
  return(solnpResult)
}
forDeltaPhaiOptim <- function(params, freq, output=FALSE) {
  delta <- params[[1]]
  phai  <- params[[2]]
  if (output == TRUE)
    solnpResult <- GNiMSModel(delta, phai, freq)
  if (output == FALSE) {
    sink(nullfile())
    solnpResult <- GNiMSModel(delta, phai, freq)
    sink()
  }
  optimizedFuncValue <- solnpResult$value[length(solnpResult$value)]
  return(optimizedFuncValue)
}
DisplayResult <- function(freq) {
  optimValues <- optim(c(1,1), forDeltaPhaiOptim , freq = freq)
  optimDelta <- optimValues$par[[1]]
  optimPhai <- optimValues$par[[2]]
  result <- GNiMSModel(optimDelta, optimPhai, freq)
  result$modelParams <- optimValues$par
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

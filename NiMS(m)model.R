rm(list = ls(all.names = TRUE))
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
rows <- CountRow(freq)
### m=1 -> count(delta)=1; m=2 -> count(delta)=2 
NiMSmConstrFunc <- function(p, m, deltas, ...) {
  rows <- CountRow(p)
  NiMSm_0Sum_Constr <- c(sum(p) - 1)
  for (i in 1:rows) {
    W1 <- sum(CalcW1Area(p, i))
    W2 <- sum(CalcW2Area(p, i))
    delta_m_i <- 1
    for (k in 1:m){
      delta_m_i <- delta_m_i * deltas[[k]]^(i^(k-1))
    }
    NiMSm_0Sum_Constr <- append(NiMSm_0Sum_Constr, W1 - delta_m_i * W2)
  }
  return (NiMSm_0Sum_Constr)
}
NiMSmModel <- function(m, deltas, freq) {
  p0 <- rep(1/length(freq), length(freq))
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  eqB <- rep(0, rows + 1)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = NiMSmConstrFunc, eqB = eqB,
                       LB = lowerBound, freq = freq, m = m, deltas = deltas)
  solnpResult$df <- rows - m
  return(solnpResult)
}
forDeltasOptim <- function(m, deltas, freq, output=FALSE) {
  if (output == TRUE)
    solnpResult <- NiMSmModel(m, deltas, freq)
  if (output == FALSE) {
    sink(nullfile())
    solnpResult <- NiMSmModel(m, deltas, freq)
    sink()
  }
  optimizedFuncValue <- solnpResult$value[length(solnpResult$value)]
  return(optimizedFuncValue)
}
DisplayNiMSmResult <- function(m) {
  deltas <- rep(1, m)
  if (m == 1) {
    optimValues <- optimize(forDeltasOptim, interval = c(0, 10), freq = freq, m = m)
    modelParams <- optimValues$minimum 
  }else{
    optimValues <- optim(deltas, forDeltasOptim , freq = freq, m = m)
    modelParams <- optimValues$par
  }
  result <- NiMSmModel(m, modelParams, freq)
  result$modelParams <- modelParams
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
m <- 1
system.time(result <- DisplayNiMSmResult(1))
result$modelParams
mhat <- result$pars * sum(freq)
mhat


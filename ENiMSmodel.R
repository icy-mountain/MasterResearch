rm(list = ls(all.names = TRUE))
source("./utilities.R")
hoge <- c(2, 3, 2, 1, #i=1
          1, 8, 6, 5, #i=2
          2, 4, 17, 5,#i=3
          1, 3, 3, 10)#i=4
freq <- hoge
## "..." is solution of nsolnp bug
ENiMSConstrFunc <- function(p, ...) {
  rows <- CountRow(p)
  ENiMS_0Sum_Constr <- c()
  ENiMS_0Sum_Constr <- append(ENiMS_0Sum_Constr, sum(p) - 1)
  for (i in 1:(rows-1)) {
    W1 <- sum(CalcW1Area(p, i))
    W2 <- sum(CalcW2Area(p, i))
    ENiMS_0Sum_Constr <- append(ENiMS_0Sum_Constr, W1 / (W2 + 1e-15))
  }
  return (ENiMS_0Sum_Constr)
}
ENiMSModel <- function(delta, freq) {
  p0 <- rep(1/length(freq), length(freq))
  paramLowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  eqB <- c(0)
  eqB <- append(eqB, rep(delta, rows-1))
  ENiMSSolnp <- solnp(p0, fun = objectFunc, eqfun = ENiMSConstrFunc, eqB = eqB,
                        LB = paramLowerBound, freq = freq)
  ENiMSSolnp$df <- rows - 1
  return(ENiMSSolnp)
}
ENiMSDeltaOptim <- function(delta, freq, output=FALSE) {
  if (output == TRUE)
    solnpResult <- ENiMSModel(delta, freq)
  if (output == FALSE) {
    sink(nullfile())
    solnpResult <- ENiMSModel(delta, freq)
    sink()
  }
  optimValue <- solnpResult$value[length(solnpResult$value)]
  return(optimValue)
}
optimValues <- optimize(ENiMSDeltaOptim, interval = c(0, 10), freq = freq)
optimDelta <- optimValues$minimum
result <- ENiMSModel(delta = optimDelta,freq = freq)
mhat <- result$pars * sum(freq)
G2 <- CalcG2(freq, mhat)
G2
pchisq(q=G2, df=result$df, lower.tail=FALSE)
CalcAICplus(freq, result)

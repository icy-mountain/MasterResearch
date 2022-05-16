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
freq <- tab3_freq
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
system.time(optimValues <- optim(c(1,1), forDeltaPhaiOptim , freq = freq))
optimDelta <- optimValues$par[[1]]
optimPhai <- optimValues$par[[2]]
optimDelta
optimPhai
optimDelta * optimPhai^2
result <- GNiMSModel(optimDelta, optimPhai, freq)
mhat <- result$pars * sum(freq)
mhat
G2 <- CalcG2(freq, mhat)
G2
pchisq(q=G2, df=result$df, lower.tail=FALSE)
CalcAICplus(freq, result)

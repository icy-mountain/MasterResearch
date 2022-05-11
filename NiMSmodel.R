rm(list = ls(all.names = TRUE))
source("./utilities.R")
hoge <- c(2, 3, 2, 1, #i=1
          1, 8, 6, 5, #i=2
          2, 4, 17, 5,#i=3
          1, 3, 3, 10)#i=4
freq <- hoge
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
result <- NiMSModel(freq = freq)
mhat <- result$pars * sum(freq)
mhat
G2 <- CalcG2(freq, mhat)
G2
pchisq(q=G2, df=result$df, lower.tail=FALSE)
CalcAICplus(freq, result)

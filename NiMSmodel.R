rm(list = ls(all.names = TRUE))
source("./utilities.R")
hoge <- c(2, 3, 2, 1, #i=1
          1, 8, 6, 5, #i=2
          2, 4, 17, 5,#i=3
          1, 3, 3, 10)#i=4
freq <- hoge
NiMSConstrFunc <- function(p, ...) {
  rows <- CountRow(p)
  NiMS_0Sum_Constr <- c()
  NiMS_0Sum_Constr <- append(NiMS_0Sum_Constr, sum(p) - 1)
  for (i in 1:(rows-1)) {
    W1 <- sum(CalcW1Area(p, i))
    W2 <- sum(CalcW2Area(p, i))
    NiMS_0Sum_Constr <- append(NiMS_0Sum_Constr, W1 - W2)
  }
  return (NiMS_0Sum_Constr)
}
NiMSModel <- function(freq) {
  p0 <- rep(1/length(freq), length(freq))
  paramLowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  NiMSSolnp <- solnp(p0, fun = objectFunc, eqfun = NiMSConstrFunc,
                     eqB = rep(0, rows),
                     LB = paramLowerBound, freq = freq)
  NiMSSolnp$df <- rows
  return(NiMSSolnp)
}
result <- NiMSModel(freq = freq)
mhat <- result$pars * sum(freq)
G2 <- CalcG2(freq, mhat)
pchisq(q=G2, df=result$df, lower.tail=FALSE)
CalcAICplus(freq, result)

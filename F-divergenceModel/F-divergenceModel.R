rm(list = ls(all.names = TRUE))
#install.packages('Rsolnp')
library(Rsolnp)
source("./utilities.R")
tab1_freq <-
  c(1520, 266, 124, 66,
    234,  1512,432, 78,
    117,  362, 1772,205,
    36,   82,  179, 492)
tab2_freq <-
  c(50, 45, 8, 18,  8,
    28, 174,84,154, 55,
    11, 78, 110,223,96,
    14, 150,185,714, 447,
    3,  42, 72, 320, 411)
# solnp section####
ASkfConstrFunc <- function(p, ft, name, score, k, ...) {
  rows <- CountRow(p)
  ASkf_0Sum_Constr <- c(sum(p) - 1)
  F2pc <- CalcF2pc(p, ft, name)
  alphas <- CalcAllAlphas(F2pc, score, k)
  if (k < r-1) {
    for (j in (k+2):rows) {
      alpha_score_sum <- CalcAlphaScoreSum(alphas, score, k, 1, j)
      IdxIJ <- RowColToIdx(p, 1, j)
      IdxJI <- RowColToIdx(p, j, 1)
      ASkf_0Sum_Constr <- append(ASkf_0Sum_Constr,
                                 F2pc[[IdxIJ]] - F2pc[[IdxJI]] - alpha_score_sum)
    }
  }
  for (i in 2:(rows-1)) {
    for (j in (i+1):rows) {
      alpha_score_sum <- CalcAlphaScoreSum(alphas, score, k, i, j)
      IdxIJ <- RowColToIdx(p, i, j)
      IdxJI <- RowColToIdx(p, j, i)
      ASkf_0Sum_Constr <- append(ASkf_0Sum_Constr,
                                 F2pc[[IdxIJ]] - F2pc[[IdxJI]] - alpha_score_sum)
    }
  }
  return (ASkf_0Sum_Constr)
}
ASkfModel <- function(freq, ft, name, score, k) {
  p0 <- rep(1/length(freq), length(freq))
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  constr_num <- rows * (rows - 1) / 2 - k
  eqB <- rep(0, constr_num + 1)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = ASkfConstrFunc, eqB = eqB,
                       LB = lowerBound, freq = freq, ft = ft, name = name,
                       score = score, k = k)
  solnpResult$df <- rows * (rows-1) / 2 - k
  return(solnpResult)
}
DisplayASkfResult <- function(freq, ft, name, score, k) {
  result <- ASkfModel(freq, ft, name, score, k)
  F2pc <- CalcF2pc(result$pars, ft, name)
  alphas <- CalcAllAlphas(F2pc, score, k)
  pHat <- result$pars
  mhat <- result$pars * sum(freq)
  result$G2 <- CalcG2(freq, mhat)
  result$pValue <- pchisq(q=result$G2, df=result$df, lower.tail=FALSE)
  result$alphas <- alphas
  print(sprintf("k:%s", k))
  print(sprintf("f:%s", format(ft)))
  print(sprintf("df:%s", result$df))
  print(sprintf("G2:%s", result$G2))
  print(sprintf("pValue:%s", result$pValue))
  for (i in 1:k) {
    print(sprintf("alpha_%d:%f", i, alphas[[i]]))
  }
  return(result)
}
# Operation check section####
freq <- tab1_freq
r <- CountRow(freq)
ft <- ~ (1 - t)^2
name <- "t"
score <- 1:r
k <- r-1
result <- DisplayASkfResult(freq, ft, name, score, k)
result$pars * sum(freq)
result$alphas

# Variance check section####
hess <- result$hessian 
invHess <- solve(hess)
F2Sigmac <- CalcF2pc(invHess, ft, name)
alphasSigma <- CalcAllAlphas(F2Sigmac, score, k)
alphasSigma

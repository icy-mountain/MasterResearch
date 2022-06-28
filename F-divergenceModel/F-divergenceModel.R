rm(list = ls(all.names = TRUE))
install.packages('Rsolnp')
library(Rsolnp)
source("./utilities.R")
# hoge ####
tab1_freq <-
  c(1520, 266, 124, 66,
    234,  1512,432, 78,
    117,  362, 1772,205,
    36,   82,  179, 492)
tab1_p <- tab1_freq / sum(tab1_freq)
# df: r * (r-1) / 2 - k
r <- CountRow(tab1_freq)
f <- ~ (1 - t)^2
name <- "t"
score <- 1:r
Calc2_pij_c(tab1_freq, 1, 3)
F2pc <- CalcF_2_p_c(tab1_p, f, name)
F2pc
alphas <- rep(0, r-1)
# i=1 j=2 score=1234 
idxIJ <- RowCol2Idx(F2pc, 1, 2)
idxJI <- RowCol2Idx(F2pc, 2, 1)
alphas[[1]] <- (F2pc[[idxIJ]] + F2pc[[idxJI]] + CalcNumerator(score, alphas, 2)) /
            CalcDenominator(score, 2)
alphas
CalcAlpha <- function(F2pc, i, j, score, alphas) {
  idxIJ <- RowCol2Idx(F2pc, i, j)
  idxJI <- RowCol2Idx(F2pc, j, i)
  num <- CalcNumerator(score, alphas, j)
  den <- CalcDenominator(score, j)
  ans <- (F2pc[[idxIJ]] + F2pc[[idxJI]] + num) / den
  print(sprintf("ij:%s ji:%s", F2pc[[idxIJ]], F2pc[[idxJI]]))
  print(sprintf("num:%s  den:%s", num, den))
  return(ans)
}
alphas[[3]] <-  CalcAlpha(F2pc, 1, 4, score, alphas)

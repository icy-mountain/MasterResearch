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
# df: r * (r-1) / 2 - k
r <- CountRow(tab1_freq)
f <- ~ (1 - t)^2
name <- "t"
score <- 1:r
Calc2_pij_c(tab1_freq, 1, 3)
p0 <- rep(exp(1), length(tab1_freq))
p0 <- tab1_freq / sum(tab1_freq)
hoge <- CalcF_2_pi_c(p0, ~ t*log(t), "t")
PrintRowsVec(hoge)
PrintDiagonal(hoge)

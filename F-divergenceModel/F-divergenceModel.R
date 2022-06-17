rm(list = ls(all.names = TRUE))
source("./utilities.R")
# hoge ####
tab1_freq <-
  c(2, 3, 2, 1,
    1, 8, 6, 5,
    2, 4, 17, 5,
    1, 3, 3, 10)
tab1_freq_T <- 
  c(2, 1, 2, 1,
    3, 8, 4, 3,
    2, 6, 17,3,
    1, 5, 5, 10)
p0 <- rep(exp(1), length(tab1_freq))
p0 <- tab1_freq / sum(tab1_freq)
hoge <- CalcF_2_pi_c(p0, ~ t*log(t), "t")
PrintRowsVec(hoge)
PrintDiagonal(hoge)

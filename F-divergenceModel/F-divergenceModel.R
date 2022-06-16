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
TransoseVec <- function(vec) {
  r <- CountRow(vec)
  transposed <- rep(1, length(vec))
  for (i in 1:r) {
    start <- (i-1)*r + 1
    end <- start + r - 1
    iRowVec <- vec[start:end]
    for (j in 1:r) {
      idx <- (j-1)*r + i
      transposed[[idx]] <- iRowVec[[j]]
    }
  }
  return(transposed)
}
printRowsVec <- function(vec) {
  r <- CountRow(vec)
  for (i in 1:r) {
    start <- (i-1)*r + 1
    end <- start + r - 1
    print(vec[start:end])
  }
}
printDiagonal <- function(vec) {
  r <- CountRow(vec)
  for(i in 1:r){
    idx <- (i-1)*r + i
    print(vec[[idx]])
  }
}
Calc2_pi_c <- function(p) {
  r <- CountRow(p)
  p_T <- TransoseVec(p)
  pc <- 2 * p / (p + p_T)
  return(pc)
}
CalcF_2_pi_c <- function(p, f_formula, name) {
  #pi_c <- Calc2_pi_c(p)
  pi_c <- p
  deriv_func <- deriv(f_formula, c(name), func = TRUE)
  evalated <- deriv_func(pi_c)
  F_2_pi_c <- attributes(evalated)$gradient
  return(F_2_pi_c)
}
p0 <- tab1_freq / sum(tab1_freq)
p0 <- rep(1, length(tab1_freq))
hoge <- CalcF_2_pi_c(p0, ~ t*log(t), "t")
hoge

CountRow <- function(freq) {
  r <- ifelse(floor(sqrt(length(freq))) < ceiling(sqrt(length(freq))), -1, sqrt(length(freq)))
  return(r)
}
Idx2RowCol <- function(idx, rows) {
  return(c((idx-1) %/% rows+1, (idx-1) %% rows+1))
}
removeZero <- function(freq) {
  return(freq[freq > 0])
}
fullModel <- function(freq) {
  nonZeroFreq <- removeZero(freq)
  nonZeroMean <- removeZero(freq / sum(freq))
  return(-sum(nonZeroFreq * log(nonZeroMean)))
}
objectFunc <- function(p, freq, ...) {
  return(-sum(freq * log(p)))
}
PrintRowsVec <- function(vec) {
  r <- CountRow(vec)
  for (i in 1:r) {
    start <- (i-1)*r + 1
    end <- start + r - 1
    print(vec[start:end])
  }
}
PrintDiagonal <- function(vec) {
  r <- CountRow(vec)
  for(i in 1:r){
    idx <- (i-1)*r + i
    print(vec[[idx]])
  }
}
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
Calc2_pi_c <- function(p) {
  p_T <- TransoseVec(p)
  pc <- 2 * p / (p + p_T)
  return(pc)
}
CalcF_2_pi_c <- function(p, f_formula, name) {
  pi_c <- Calc2_pi_c(p)
  deriv_func <- deriv(f_formula, c(name), func = TRUE)
  evalated <- deriv_func(pi_c)
  F_2_pi_c <- attributes(evalated)$gradient
  return(F_2_pi_c)
}
CalcG2 <- function(freq, mhat) {
  return(2 * sum(freq * (log(freq) - log(mhat))))
}
CalcAICplus <- function(freq, solnpResult) {
  mhat <- solnpResult$pars * sum(freq)
  G2 <- CalcG2(freq, mhat)
  return(G2 - 2 * solnpResult$df)
}
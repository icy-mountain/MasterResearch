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
# W1 = P(X<=i, Y>=i) 
ExtractW1Area <- function(p, i) {
  IdxVec <- 1:length(p)
  W1Indices <- c()
  for (idx in IdxVec) {
    r_c <- Idx2RowCol(idx, CountRow(p))
    if (r_c[[1]] <= i && r_c[[2]] >= i)
      W1Indices <- append(W1Indices, idx)
  }
  return(p[W1Indices])
}
# W2 = P(X>=i, Y<=i) 
ExtractW2Area <- function(p, i) {
  IdxVec <- 1:length(p)
  W1Indices <- c()
  for (idx in IdxVec) {
    r_c <- Idx2RowCol(idx, CountRow(p))
    if (r_c[[1]] >= i && r_c[[2]] <= i)
      W1Indices <- append(W1Indices, idx)
  }
  return(p[W1Indices])
}
CalcG2 <- function(freq, mhat) {
  return(2 * sum(freq * (log(freq) - log(mhat))))
}
CalcAICplus <- function(freq, solnpResult) {
  mhat <- solnpResult$pars * sum(freq)
  G2 <- CalcG2(freq, mhat)
  return(G2 - 2 * solnpResult$df)
}
objectFunc <- function(p, freq, ...) {
  return(-sum(freq * log(p)))
}
SumAllW1Area <- function(p) {
  ans <- 0
  rows <- CountRow(p)
  for (i in 1:rows)
    ans <- ans + sum(ExtractW1Area(p, i))
  return(ans)
}
SumAllW2Area <- function(p) {
  ans <- 0
  rows <- CountRow(p)
  for (i in 1:rows)
    ans <- ans + sum(ExtractW2Area(p, i))
  return(ans)
}
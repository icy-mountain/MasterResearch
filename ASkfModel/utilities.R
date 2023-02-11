CountRow <- function(freq) {
  r <- ifelse(floor(sqrt(length(freq))) < ceiling(sqrt(length(freq))), -1, sqrt(length(freq)))
  return(r)
}
IdxToRowCol <- function(idx, rows) {
  return(c((idx-1) %/% rows+1, (idx-1) %% rows+1))
}
RowColToIdx <- function(vec, row, col) {
  r <- CountRow(vec)
  return((row - 1) * r + col)
}
P_ <- function(p, i, j){
  return(p[[RowColToIdx(p, i, j)]])
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
CalcG2 <- function(freq, mhat) {
  return(2 * sum(freq * (log(freq) - log(mhat))))
}
CalcAICplus <- function(freq, solnpResult) {
  mhat <- solnpResult$pars * sum(freq)
  G2 <- CalcG2(freq, mhat)
  return(G2 - 2 * solnpResult$df)
}
CalcVandermondeMatrix <- function(score, k) {
  pow <- function(a, n) {
    return(a^n)
  }
  ans <- matrix(1, nrow=1, ncol=k)
  if (k == 1)
    return(ans)
  for (i in 1:(k-1)) {
    ans <- rbind(ans, unlist(map(score[1:k], pow, i)))
  }
  return(ans)
}
CalcVandermondeMatrix2 <- function(score, k) {
  if (k == 1)
    return(1)
  pow <- function(a, n) {
    return(a^n)
  }
  ans <- matrix(score[2:(k+1)] - rep(score[[1]], k), nrow=1, ncol=k)
  for (i in 2:k) {
    ans <- rbind(ans, unlist(map(score[2:(k+1)], pow, i)) - unlist(map(rep(score[[1]], k), pow, i)))
  }
  return(ans)
}
CalcAntiVandermondeMatrix <- function(score, k) {
  for (i in 1:k) {
    rowVec <- matrix(0, nrow=1, ncol=k)
    for (j in 1:i) {
      value <- 1
      if (j > 1) {
        for (l in 1:(j-1)) {
          # value <- value * score[[i - l]] #ここが難しい
          # value <- value * (score[[i]] - score[[l]])
          value <- value * (i - l)
        }
      }
      # rowVec[[j]] <- score[[i]] * value
      rowVec[[j]] <- i * value
    }
    if (i == 1)
      ans <- rowVec
    else
      ans <- rbind(ans, rowVec)
  }
  return(ans)
}

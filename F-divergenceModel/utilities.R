CountRow <- function(freq) {
  r <- ifelse(floor(sqrt(length(freq))) < ceiling(sqrt(length(freq))), -1, sqrt(length(freq)))
  return(r)
}
Idx2RowCol <- function(idx, rows) {
  return(c((idx-1) %/% rows+1, (idx-1) %% rows+1))
}
RowCol2Idx <- function(vec, row, col) {
  r <- CountRow(vec)
  return((row - 1) * r + col)
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
Calc2_pij_c <- function(p, i, j) {
  r <- CountRow(p)
  pij <- p[[RowCol2Idx(p, i, j)]]
  pji <- p[[RowCol2Idx(p, j, i)]]
  return(2 * pij / (pij + pji))
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
Calc2_p_c <- function(p) {
  p_T <- TransoseVec(p)
  pc <- 2 * p / (p + p_T)
  return(pc)
}
CalcF_2_p_c <- function(p, f_formula, name) {
  pi_c <- Calc2_p_c(p)
  deriv_func <- deriv(f_formula, c(name), func = TRUE)
  evalated <- deriv_func(pi_c)
  F_2_pi_c <- attributes(evalated)$gradient
  return(F_2_pi_c)
}
CalcDenominator <- function(score, n) {
  ans <- 1
  for (i in 1:(n-1)) {
    ans <- ans * (score[[n]] - score[[i]])
  }
  return(ans)
}
CalcNumerator <- function(score, alpha, n) {
  ans <- 0
  i <- 1
  while(i <= (n-2)) {
    #print(sprintf("i:%s n:%s", i, n))
    #print(sprintf("called!"))
    mul <- 1
    for (j in 1:i) {
      mul <- mul * (score[[n]] - score[[j]])
    }
    ans <- ans + mul * alpha[[i]]
    i <- i + 1
  }
  return(ans)
}
CalcAlpha <- function(F2pc, i, j, score, alphas) {
  idxIJ <- RowCol2Idx(F2pc, i, j)
  idxJI <- RowCol2Idx(F2pc, j, i)
  num <- CalcNumerator(score, alphas, j)
  den <- CalcDenominator(score, j)
  ans <- (F2pc[[idxIJ]] - F2pc[[idxJI]] + num) / (-1 * den)
  return(ans)
}
CalcAllAlphas <- function(F2pc, score, k) {
  alphas <- rep(0, k)
  for (i in 1:k) {
    alphas[[i]] <- CalcAlpha(F2pc, 1, i+1, score, alphas)
  }
  return(alphas)
}
CalcAlphaScoreSum <- function(alphas, score, k, i, j) {
  alpha_score_sum <- 0
  for (h in 1:k) {
    score_sum <- 0
    score_mul <- 1
    # left
    for (g in 1:h) {
      score_mul <- score_mul * (score[[i]] - score[[g]])
    }
    score_sum <- score_sum + score_mul
    score_mul <- 1
    # right
    for (g in 1:h) {
      score_mul <- score_mul * (score[[j]] - score[[g]])
    }
    score_sum <- score_sum - score_mul
    alpha_score_sum <- alpha_score_sum + score_sum * alphas[[h]]
  }
  return(alpha_score_sum)
}
CalcG2 <- function(freq, mhat) {
  return(2 * sum(freq * (log(freq) - log(mhat))))
}
CalcAICplus <- function(freq, solnpResult) {
  mhat <- solnpResult$pars * sum(freq)
  G2 <- CalcG2(freq, mhat)
  return(G2 - 2 * solnpResult$df)
}
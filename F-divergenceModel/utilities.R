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
Calc2pijc <- function(p, i, j) {
  r <- CountRow(p)
  pij <- p[[RowColToIdx(p, i, j)]]
  pji <- p[[RowColToIdx(p, j, i)]]
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
Calc2pc <- function(p) {
  p_T <- TransoseVec(p)
  pc <- 2 * p / (p + p_T)
  return(pc)
}
CalcF2pc <- function(p, f_formula, name) {
  pi_c <- Calc2pc(p)
  deriv_func <- deriv(f_formula, c(name), func = TRUE)
  evalated <- deriv_func(pi_c)
  F_2_pi_c <- attributes(evalated)$gradient
  return(F_2_pi_c)
}
CalcDerivF2pc <- function(p, f_formula, name) {
  pi_c <- Calc2pc(p)
  deriv_func <- deriv(f_formula, c(name), func = TRUE, hessian = TRUE)
  evalated <- deriv_func(pi_c)
  derivF_2_pi_c <- attributes(evalated)$hessian
  return(derivF_2_pi_c)
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
  idxIJ <- RowColToIdx(F2pc, i, j)
  idxJI <- RowColToIdx(F2pc, j, i)
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
    multiply_left <- 1
    multiply_right <- 1
    for (g in 1:h) {
      multiply_left <- multiply_left * (score[[i]] - score[[g]])
      multiply_right <- multiply_right * (score[[j]] - score[[g]])
    }
    score_sum <- multiply_left - multiply_right
    alpha_score_sum <- alpha_score_sum + score_sum * alphas[[h]]
  }
  return(alpha_score_sum)
}
CalcDerivAlpha <- function(funcF, p, i, j, score, derivAlphas) {
  pij <- p[[RowColToIdx(p, i, j)]]
  pji <- p[[RowColToIdx(p, j, i)]]
  num_dPij <- CalcNumerator(score, derivAlphas[["dPij"]], j)
  num_dPji <- CalcNumerator(score, derivAlphas[["dPji"]], j)
  den <- CalcDenominator(score, j)
  exprAlpha_dPij <- (funcF(2 * pij / (pij + pji)) -
                       funcF(2 * pji / (pij + pji)) +
                       num_dPij) /
    (-1 * den) ~ pij
  exprAlpha_dPji <- (funcF(2 * pij / (pij + pji)) -
                       funcF(2 * pji / (pij + pji)) +
                       num_dPji) /
    (-1 * den) ~ pji
  alpha_dPij <- D(exprAlpha_dPij, pji = pji, num_dPij = num_dPij, den = den)
  alpha_dPji <- D(exprAlpha_dPji, pij = pij, num_dPji = num_dPji, den = den)
  derivAlphas[["dPij"]] <- append(derivAlphas[["dPij"]], alpha_dPij(pij = pij))
  derivAlphas[["dPji"]] <- append(derivAlphas[["dPji"]], alpha_dPji(pji = pji))
  return(derivAlphas)
}
CalcDerivAllAlphas <- function(funcF, p, score, k) {
  derivAlphas <- vector("list", length = 2)
  names(derivAlphas) <- c("dPij","dPji")
  for (i in 1:k) {
    derivAlphas <- CalcDerivAlpha(funcF, p, 1, i+1, score, derivAlphas)
  }
  return(derivAlphas)
}
MakeAlpha_dPMatrix <- function(alphas_dP, r, k){
  row <- r * (r - 1) / 2 - k + 1
  col <- r * r
  size <- length(alphas_dP[["dPij"]])
  ans <- rep(0, row * col)
  dim(ans) <- c(row, col)
  ans[1,][2:(size+1)] <- alphas_dP[["dPij"]]
  for (idx in 2:row){
    ans[idx,][1] <- alphas_dP[["dPji"]][[idx-1]]
  }
  return(ans)
}
CalcG2 <- function(freq, mhat) {
  return(2 * sum(freq * (log(freq) - log(mhat))))
}
CalcAICplus <- function(freq, solnpResult) {
  mhat <- solnpResult$pars * sum(freq)
  G2 <- CalcG2(freq, mhat)
  return(G2 - 2 * solnpResult$df)
}
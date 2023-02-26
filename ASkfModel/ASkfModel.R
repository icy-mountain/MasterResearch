# Note: ####
# alpha in signature means alpha_star.
# Sigma -> variance-covariance matrix of pi.
# AlphaSigma -> (dalpha/dpi) %*% Sigma %*% (dalpha/dpi)^T
# solnp section ############
Calc2pc <- function(p) {
  p_T <- TransoseVec(p)
  pc <- 2 * p / (p + p_T)
  return(pc)
}
CalcF2pc <- function(p, f, name) {
  pi_c <- Calc2pc(p) + (.Machine$double.eps)^(1/3) #1.0e-6 # for log(non-zero)
  f_formula <- as.formula(paste(name, " ~ " ,f))
  deriv_func <- deriv(f_formula, c(name), func = TRUE)
  evaluated <- deriv_func(pi_c)
  F_2_pi_c <- attributes(evaluated)$gradient
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
  # if for S model
  if (k < 1)
    return(0)
  alphas <- rep(0, k)
  for (i in 1:k) {
    alphas[[i]] <- CalcAlpha(F2pc, 1, i+1, score, alphas)
  }
  return(alphas)
}
CalcAlphaScoreSum <- function(alphas, score, k, i, j) {
  # if for S model
  if (k < 1)
    return(0)
  alpha_score_sum <- 0
  for (h in 1:k) {
    multiply_left <- 1
    multiply_right <- 1
    for (g in 1:h) {
      multiply_left <- multiply_left * (score[[i]] - score[[g]])
      multiply_right <- multiply_right * (score[[j]] - score[[g]])
    }
    score_sum <- multiply_left - multiply_right
    alpha_score_sum <- alpha_score_sum + (score_sum * alphas[[h]])
  }
  return(alpha_score_sum)
}
ASkfConstrFunc <- function(p, ASkf_f, name, score, k, ...) {
  rows <- CountRow(p)
  ASkf_0Sum_Constr <- c(sum(p) - 1)
  F2pc <- CalcF2pc(p, ASkf_f, name)
  alphas <- CalcAllAlphas(F2pc, score, k)
  if (k < rows-1) {
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
ASkfModel <- function(freq, f, name, score, k, ctrl) {
  p0 <- rep(1/length(freq), length(freq))
  lowerBound <- rep(0, length(freq))
  rows <- CountRow(freq)
  constr_num <- rows * (rows - 1) / 2 - k
  eqB <- rep(0, constr_num + 1)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = ASkfConstrFunc, eqB = eqB,
                       LB = lowerBound,control = ctrl, freq = freq, ASkf_f = f, name = name,
                       score = score, k = k)
  solnpResult$df <- rows * (rows - 1) / 2 - k
  return(solnpResult)
}
# covarience section ############
CalcAlphaSigma <- function(freq, f, name, score, k, ASkf_result) {
  params <- MakeASkfParamList(f, name, score, k)
  params$freq <- freq
  params$rows <- CountRow(freq)
  ASkf_Sigma <- CalcSigma(params, ASkf_result)
  f_formula <- as.formula(paste(params$f, "~ ", name))
  funcF <- makeFun(D(f_formula))
  p <- ASkf_result$pars
  alpha_dPMat <- MakeAlphaDiffPMatrix(f, name, p, params$rows, k, score)
  alphaSigma <- alpha_dPMat %*% ASkf_Sigma %*% t(alpha_dPMat)
  return(alphaSigma)
}
CalcSigma <- function(params, ASkf_result){
  p <- ASkf_result$pars 
  Dp <- diag(p)
  H <- NumDiff(ASkfConstr_No0Sum, p, params)
  inv_HDH <- solve(t(H) %*% Dp %*% H)
  Sigma <- Dp - p %*% t(p) - Dp %*% H %*% inv_HDH %*% t(H) %*% Dp
  return(Sigma)
}
NumDiff <- function(hFunc, p, params) {
  rows <- CountRow(p)
  constr_num <- rows * (rows - 1) / 2 - params$k
  ans <- matrix(c(1:(constr_num * rows * rows)), nrow = rows * rows)
  eps <- (.Machine$double.eps)^(1/3)
  d <- eps * p + eps  
  E <- diag(d)
  for (j in 1:length(p)) {
    f1 <- hFunc(p + E[,j], params)
    f2 <- hFunc(p - E[,j], params)
    Ft <- (f1 - f2) / (2 * d[j])
    ans[j,] <- Ft
  }
  return(ans)
}
# partial diff of alpha by pi
PartDiff <- function(f, name, p, i, j, k, score){
  minimum <- (.Machine$double.eps)^(1/3)
  p_plus <- p
  p_minus <- p
  p_plus[[RowColToIdx(p_plus, i, j)]] <- P_(p_plus, i, j) + minimum
  p_minus[[RowColToIdx(p_minus, i, j)]] <- P_(p_minus, i, j) - minimum
  F2pc_plus <- CalcF2pc(p_plus, f, name)
  F2pc_minus <- CalcF2pc(p_minus, f, name)
  return((CalcAllAlphas(F2pc_plus, score, k)
          - CalcAllAlphas(F2pc_minus, score, k))
         / (2 * minimum))
}
# make matrix by PartDiff
MakeAlphaDiffPMatrix <- function(f, name, p, r, k, score){
  row <- k
  col <- r * r
  ans <- rep(0, row * col)
  dim(ans) <- c(row, col)
  for (colIdx in 1:col){
    ij <- IdxToRowCol(colIdx, r)
    i <- ij[[1]]
    j <- ij[[2]]
    ans[,colIdx] <- PartDiff(f, name, p, i, j, k, score)
  }
  return(ans)
}
ASkfConstr_No0Sum <- function(p, params, ...) {
  rows <- CountRow(p)
  f <- params$f
  name <- params$name
  score <- params$score
  k <- params$k
  ASkf_Constr <- c()
  F2pc <- CalcF2pc(p, f, name)
  alphas <- CalcAllAlphas(F2pc, score, k)
  if (k < rows-1) {
    for (j in (k+2):rows) {
      alpha_score_sum <- CalcAlphaScoreSum(alphas, score, k, 1, j)
      IdxIJ <- RowColToIdx(p, 1, j)
      IdxJI <- RowColToIdx(p, j, 1)
      ASkf_Constr <- append(ASkf_Constr,
                            F2pc[[IdxIJ]] - F2pc[[IdxJI]] - alpha_score_sum)
    }
  }
  for (i in 2:(rows-1)) {
    for (j in (i+1):rows) {
      alpha_score_sum <- CalcAlphaScoreSum(alphas, score, k, i, j)
      IdxIJ <- RowColToIdx(p, i, j)
      IdxJI <- RowColToIdx(p, j, i)
      ASkf_Constr <- append(ASkf_Constr,
                            F2pc[[IdxIJ]] - F2pc[[IdxJI]] - alpha_score_sum)
    }
  }
  return (ASkf_Constr)
}
MakeASkfParamList <- function(f, name, score, k){
  params <- c()
  params$f <- f
  params$name <- name
  params$score <- score
  params$k <- k
  return(params)
}
# display section ############
FormatAlphaSigma <- function(alphas, alphaStdError, lower, upper){
  cat("alpha_stars:", "\n")
  cat("              ", "Estimate  ", "Std.Error  ", "Confidential.Interval", "\n", sep="")
  k <- length(alphas)
  for (i in 1:k) {
    cat("alpha_star", i, " ", format(c(alphas[[i]], alphaStdError[[i]]),
                                     justify = "right", digits=4, nsmall=4, width=10), sep="")
    cat("   [")
    cat(format(c(lower[[i]], upper[[i]]),
               justify = "right", digits=4, nsmall=4, width=10), sep=",")
    cat("]", sep="")
    if (lower[[i]] * upper[[i]] > 0)
      cat(" *", sep="")
    # cat(abs(upper[[i]] - lower[[i]]))
    cat("\n")
  }
  cat("(*) means interval excluding 0.", "\n")
}
FormatValues <- function(k, f, result){
  cat("**********result**********", "\n")
  cat("k:", k, "\n")
  cat("f:", format(f), "\n")
  cat("df:", result$df, "\n")
  cat("G2:", format(result$G2, digits=4, width=4), "\n")
  cat("pValue:", format(result$pValue, digits=4, width=4), "\n")
}
DisplayASkfResult <- function(freq, f, name, score, k, output = TRUE, ctrl = list(trace = 0)) {
  result <- ASkfModel(freq, f, name, score, k, ctrl)
  pHat <- result$pars
  F2pc <- CalcF2pc(pHat, f, name)
  mHat <- pHat * sum(freq)
  result$G2 <- CalcG2(freq, mHat)
  result$pValue <- pchisq(q=result$G2, df=result$df, lower.tail=FALSE)
  # if for S model
  if (k <= 0) {
    if (output) 
      FormatValues(k, f, result)
    return(result)
  }
  alphas <- CalcAllAlphas(F2pc, score, k)
  alphaStdError <- sqrt(diag(CalcAlphaSigma(freq, f, name, score, k, result)) / sum(freq))
  interval <- 1.96 * alphaStdError
  
  upper <- alphas + interval[1:k]
  lower <- alphas - interval[1:k]
  result$alphas <- alphas
  if (output) {
    FormatValues(k, f, result)
    FormatAlphaSigma(alphas, alphaStdError, lower, upper)
  }
  
  return(result)
}
# Note: ############
# p means pi:π in codes.
# alpha in signature means alpha_star.
# Sigma means variance-covariance matrix of pi.
# AlphaSigma -> (dalpha/dpi) %*% Sigma %*% (dalpha/dpi)^T
# solnp section ############
Calc2pc <- function(p) {
  p_T <- TransoseVec(p)
  pc <- 2 * p / (p + p_T)
  return(pc)
}
CalcF2pc <- function(p, f, name) {
  pi_c <- Calc2pc(p) + (.Machine$double.eps)^(1/3) # for log(non-zero)
  f_formula <- as.formula(paste(name, " ~ " ,f))
  deriv_func <- deriv(f_formula, c(name), func = TRUE)
  evaluated <- deriv_func(pi_c)
  F_2_pi_c <- attributes(evaluated)$gradient
  return(F_2_pi_c)
}
CalcDenominator <- function(score, j) {
  ans <- 1
  for (i in 1:(j-1)) {
    ans <- ans * (score[[j]] - score[[i]])
  }
  return(ans)
}
CalcNumerator <- function(score, alpha, j) {
  ans <- 0
  h <- 1
  while(h <= (j-2)) {
    mul <- 1
    for (i in 1:h) {
      mul <- mul * (score[[j]] - score[[i]])
    }
    ans <- ans + mul * alpha[[h]]
    h <- h + 1
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
  r <- CountRow(p)
  ASkf_0Sum_Constr <- c(sum(p) - 1)
  F2pc <- CalcF2pc(p, ASkf_f, name)
  alphas <- CalcAllAlphas(F2pc, score, k)
  if (k < r-1) {
    for (j in (k+2):r) {
      alpha_score_sum <- CalcAlphaScoreSum(alphas, score, k, 1, j)
      IdxIJ <- RowColToIdx(p, 1, j)
      IdxJI <- RowColToIdx(p, j, 1)
      ASkf_0Sum_Constr <- append(ASkf_0Sum_Constr,
                                 F2pc[[IdxIJ]] - F2pc[[IdxJI]] - alpha_score_sum)
    }
  }
  for (i in 2:(r-1)) {
    for (j in (i+1):r) {
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
  r <- CountRow(freq)
  constr_num <- r * (r - 1) / 2 - k
  eqB <- rep(0, constr_num + 1)
  solnpResult <- solnp(p0, fun = objectFunc, eqfun = ASkfConstrFunc, eqB = eqB,
                       LB = lowerBound,control = ctrl, freq = freq, ASkf_f = f, name = name,
                       score = score, k = k)
  solnpResult$df <- r * (r - 1) / 2 - k
  return(solnpResult)
}
# covarience section ############
CalcAlphaSigma <- function(freq, f, name, score, k, result) {
  params <- MakeASkfParamList(freq, f, name, score, k)
  Sigma <- CalcSigma(params, result)
  alpha_dPMat <- MakeAlphaDiffPMatrix(f, name, result$pars, params$r, k, score)
  alphaSigma <- alpha_dPMat %*% Sigma %*% t(alpha_dPMat)
  return(alphaSigma)
}
CalcSigma <- function(params, result){
  p <- result$pars 
  Dp <- diag(p)
  H <- NumDiff(ASkfConstr_No0Sum, p, params)
  inv_HDH <- solve(t(H) %*% Dp %*% H)
  Sigma <- Dp - p %*% t(p) - Dp %*% H %*% inv_HDH %*% t(H) %*% Dp
  return(Sigma)
}
NumDiff <- function(hFunc, p, params) {
  r <- CountRow(p)
  constr_num <- r * (r - 1) / 2 - params$k
  ans <- matrix(c(1:(constr_num * r * r)), nrow = r * r)
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
  eps <- (.Machine$double.eps)^(1/3)
  p_plus <- p
  p_minus <- p
  p_plus[[RowColToIdx(p_plus, i, j)]] <- P_(p_plus, i, j) + eps
  p_minus[[RowColToIdx(p_minus, i, j)]] <- P_(p_minus, i, j) - eps
  F2pc_plus <- CalcF2pc(p_plus, f, name)
  F2pc_minus <- CalcF2pc(p_minus, f, name)
  return((CalcAllAlphas(F2pc_plus, score, k)
          - CalcAllAlphas(F2pc_minus, score, k))
         / (2 * eps))
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
  r <- CountRow(p)
  f <- params$f
  name <- params$name
  score <- params$score
  k <- params$k
  ASkf_Constr <- c()
  F2pc <- CalcF2pc(p, f, name)
  alphas <- CalcAllAlphas(F2pc, score, k)
  if (k < r-1) {
    for (j in (k+2):r) {
      alpha_score_sum <- CalcAlphaScoreSum(alphas, score, k, 1, j)
      IdxIJ <- RowColToIdx(p, 1, j)
      IdxJI <- RowColToIdx(p, j, 1)
      ASkf_Constr <- append(ASkf_Constr,
                            F2pc[[IdxIJ]] - F2pc[[IdxJI]] - alpha_score_sum)
    }
  }
  for (i in 2:(r-1)) {
    for (j in (i+1):r) {
      alpha_score_sum <- CalcAlphaScoreSum(alphas, score, k, i, j)
      IdxIJ <- RowColToIdx(p, i, j)
      IdxJI <- RowColToIdx(p, j, i)
      ASkf_Constr <- append(ASkf_Constr,
                            F2pc[[IdxIJ]] - F2pc[[IdxJI]] - alpha_score_sum)
    }
  }
  return (ASkf_Constr)
}
MakeASkfParamList <- function(freq, f, name, score, k){
  params <- c()
  params$freq <- freq
  params$f <- f
  params$name <- name
  params$score <- score
  params$k <- k
  params$r <- CountRow(freq)
  return(params)
}
# display section ############
FormatValues <- function(k, f, result){
  cat("**********result**********", "\n")
  cat("k:", k, "\n")
  cat("f:", format(f), "\n")
  cat("df:", result$df, "\n")
  cat("G2:", format(result$G2, digits=4, width=4), "\n")
  cat("pValue:", format(result$pValue, digits=4, width=4), "\n")
}
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
    cat("\n")
  }
  cat("(*) means interval excluding 0.", "\n")
}
FormatAlphaSigmaWithZ <- function(alphas, alphaStdError, zValue){
  cat("alpha_stars:", "\n")
  cat("              ", "Estimate  ", "Std.Error  ", "z value  ", "Pr(>|z|)", "\n", sep="")
  k <- length(alphas)
  for (i in 1:k) {
    Pr_z <- pchisq(q=zValue[[i]]^2, df=1, lower.tail=FALSE)
    if (Pr_z < 0.0001) {
      cat("alpha_star", i, " ", format(c(alphas[[i]], alphaStdError[[i]], zValue[[i]]),
                                       justify = "right", digits=4, nsmall=4, width=10), "   <0.0001", sep="")
    } else {
      cat("alpha_star", i, " ", format(c(alphas[[i]], alphaStdError[[i]], zValue[[i]], Pr_z),
                                       justify = "right", digits=4, nsmall=4, width=10), sep="")
    }
    if (0 <= Pr_z && Pr_z < 0.001)
      cat(" ***", sep="")
    if (0.001 <= Pr_z && Pr_z < 0.01)
      cat(" **", sep="")
    if (0.01 <= Pr_z && Pr_z < 0.05)
      cat(" *", sep="")
    if (0.05 <= Pr_z && Pr_z < 0.1)
      cat(" .", sep="")
    cat("\n")
  }
  cat("\n")
  cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1", "\n")
}
DisplayASkfResult <- function(freq, f, name, score, k, output = TRUE, ctrl = list(trace = 0)) {
  freq[freq <= 0] <- 0.01
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
DisplayASkfResultWithZ <- function(freq, f, name, score, k, output = TRUE, ctrl = list(trace = 0)) {
  freq[freq <= 0] <- 0.01
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
  zValue <- alphas / alphaStdError
  result$alphas <- alphas
  if (output) {
    FormatValues(k, f, result)
    FormatAlphaSigmaWithZ(alphas, alphaStdError, zValue)
  }
  
  return(result)
}
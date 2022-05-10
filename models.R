##### delete existing objects #####
rm(list = ls(all.names = TRUE))


##### data samples #####
### Yamamoto-Tomizawa2010 Table1 ###
freq <- c(29, 3, 3, 4, 5, 0, 1, 1, 9, 0, 2, 0, 7, 3, 1, 0)

### Yamamoto-Tomizawa2010 Table2 ###
freq2 <- c(6, 2, 3, 1,
           9, 4, 2, 1,
           9, 2, 3, 1,
           12, 1, 2, 1)

### Tominaga1979 (Tahata2016) ###
freq3 <- c(374, 602, 170, 64, 18, 255, 139, 71, 4, 23, 42, 55, 2, 6, 17, 53)

### Smith2006, cross-classification of GSS (Tahata-Sudo-Arimoto) ###
freq4 <- c(98, 150, 135, 53, 37, 131, 133, 43, 9, 16, 33, 15, 4, 1, 4, 21)


# set by model().
# globalAnalysResults have string ex1) "LSI1" ex2) "SU"
# r means count of rows. Used for check of square contingency table.
globalAnalysResults <- inputData <- list()
r <- 0

make_parts <- function(r) {
  array1 <- array(0, dim = c(r^2, r - 1))
  k <- 1
  for (i in 1:r) {
    for (j in 1:r) {
      if (i <= (r - 1)) array1[k, i] <- array1[k, i] + 1
      if (j <= (r - 1)) array1[k, j] <- array1[k, j] + 1
      k <- k + 1
    }
  }
  
  array2 <- c(1:r %x% 1:r)
  
  array2star <- array2
  for (i in 1:r) array2star[i + r * (i - 1)] <- 0
  
  array3 <- list()
  for (k in 1:(r - 1)) {
    tmp <- c()
    for (i in 1:r) {
      for (j in 1:r) {
        tmp <- c(tmp, j^k - i^k)
      }
    }
    array3[[k]] <- tmp
  }
  
  array4 <- array(0, dim = c(r^2, r))
  for (i in 1:r) {
    for (j in i:r) {
      array4[j + r * (j - 1), i] <- 1
      break
    }
  }
  
  array5 <- array(0, dim = c(r^2, r * (r + 1) / 2))
  k <- 1
  for (i in 1:r) {
    for (j in 1:r) {
      if (i > j) {
        array5[k, (j - 1) * (2 * r - j) / 2 + i] <- array5[k, (j - 1) * (2 * r - j) / 2 + i] + 1
      } else {
        array5[k, (i - 1) * (2 * r - i) / 2 + j] <- array5[k, (i - 1) * (2 * r - i) / 2 + j] + 1
      }
      k <- k + 1
    }
  }
  return (list(array1, array2, array3, array4, array5))
}

model <- function(freq, sort = FALSE) {
  r <<- ifelse(floor(sqrt(length(freq))) < ceiling(sqrt(length(freq))), stop(), sqrt(length(freq)))
  inputData <<- freq
  
  ##### make parts of design matrices #####
  parts <- make_parts(r)
  array1 <- parts[[1]]
  array2 <- parts[[2]]
  array3 <- parts[[3]]
  array4 <- parts[[4]]
  array5 <- parts[[5]]
  array2star <- array2
  for (i in 1:r) array2star[i + r * (i - 1)] <- 0
  
  ##### define design matrices #####
  ### SI ###
  arraySI <- array1
  ### SU ###
  arraySU <- cbind(array1, array2)
  ### SQI ###
  arraySQI <- cbind(array1, array4)
  ### SQU ###
  arraySQU <- cbind(array1, array4, array2star)
  ### S ###
  arrayS <- array5
  ### LSIk ###
  bindedArray3 <- c()
  for (i in 1:(r - 1)) {
    bindedArray3 <- cbind(bindedArray3, array3[[i]])
    assign(paste0("arrayLSI", i), cbind(array1, bindedArray3))
  }
  ### LSUk ###
  bindedArray3 <- c()
  for (i in 1:(r - 1)) {
    bindedArray3 <- cbind(bindedArray3, array3[[i]])
    assign(paste0("arrayLSU", i), cbind(array1, bindedArray3, array2))
  }
  ### LSQIk ###
  bindedArray3 <- c()
  for (i in 1:(r - 1)) {
    bindedArray3 <- cbind(bindedArray3, array3[[i]])
    assign(paste0("arrayLSQI", i), cbind(array1, bindedArray3, array4))
  }
  ### LSQUk ###
  bindedArray3 <- c()
  for (i in 1:(r - 1)) {
    bindedArray3 <- cbind(bindedArray3, array3[[i]])
    assign(paste0("arrayLSQU", i), cbind(array1, bindedArray3, array4, array2star))
  }
  ### LSk ###
  bindedArray3 <- c()
  for (i in 1:(r - 1)) {
    bindedArray3 <- cbind(bindedArray3, array3[[i]])
    assign(paste0("arrayLS", i), cbind(array5, bindedArray3))
  }

  ##### analyze with each model #####
  analysResults <- list()

  ### models without MEk & Cov=0 ###
  models <- c("SI", "SU", "SQI", "SQU", "S")
  models_k <- c("LSI", "LSU", "LSQI", "LSQU", "LS")
  for (idx in 1:length(models)) {
    formula <- as.formula(paste0("freq~array", models[idx]))
    glm <- glm(formula, family = poisson, data = list(freq))
    analysResults <- append(analysResults, list(glm))
    names(analysResults)[length(analysResults)] <- models[idx]
  }
  for (idx in 1:length(models_k)){
    for (k in 1:(r - 1)) {
      formula <- as.formula(paste0("freq~array", models_k[idx], k))
      glm <- glm(formula, family = poisson, data = list(freq))
      analysResults <- append(analysResults, list(glm))
      names(analysResults)[length(analysResults)] <- paste0(models_k[idx], k)
    }
  }

  #### MEk & Cov=0 model ####
  library(Rsolnp)
  
  #### define constraint functions ####
  ### MEk constraint###
  MEkConstrFunc <- function(p) {
    MEkConstraints <- c()
    MEkConstraints <- append(MEkConstraints, sum(p) - 1)
    
    k <- length(MEkSolnp) + 1
    for (l in 1:k) {
      momentConstraints <- c()
      for (i in 1:r) {
        for (j in 1:r) {
          momentConstraints <- append(momentConstraints, j^l - i^l)
        }
      }
      MEkConstraints <- append(MEkConstraints, sum(momentConstraints * p))
    }
    return(MEkConstraints)
  }
  ### Cov=0 constraint###
  cov0ConstrFunc <- function(p) {
    cov0Constraints <- c()
    cov0Constraints <- append(cov0Constraints, sum(p) - 1)
    
    pMatrix <- matrix(p, r, r, T)
    prodOfExp <- c()
    for (i in 1:r) {
      for (j in 1:r) {
        prodOfExp <- append(prodOfExp, rowSums(pMatrix)[i] * colSums(pMatrix)[j])
      }
    }
    cov0Constraints <- append(cov0Constraints, sum(c(1:r %x% 1:r) * (p - prodOfExp)))
    return(cov0Constraints)
  }
  #### define some functions and arguments for solnp & functions####
  removeZero <- function(freq) {
    return(freq[freq > 0])
  }
  objectFunc <- function(p) {
    return(-sum(freq * log(p)))
  }
  paramLowerBound <- rep(0, length(freq))
  equationValue <- list()
  for (i in 1:(r - 1)) equationValue <- append(equationValue, list(rep(0, i + 1)))
  p0 <- rep(1 / length(freq), length(freq))# uniform prob
  # rep(1, length.out = 5)
  # => 1 1 1 1 1

  fullModel <- function(freq) {
    nonZeroFreq <- removeZero(freq)
    nonZeroMean <- removeZero(freq / sum(freq))
    return(-sum(nonZeroFreq * log(nonZeroMean)))
  }
  moleculeOfConst <- denominatorOfConst <- 0
  for (i in 1:sum(freq)) moleculeOfConst <- moleculeOfConst + log(i)
  for (i in removeZero(freq)) {
    for (j in 1:i) {
      denominatorOfConst <- denominatorOfConst + log(j)
    }
  }
  ###MEk Solnp###
  MEkSolnp <- list()
  for (i in 1:(r - 1)) {
    tmp <- solnp(p0, fun = objectFunc, eqfun = MEkConstrFunc, eqB = equationValue[[i]], LB = paramLowerBound)
    MEkSolnp <- append(MEkSolnp, list(tmp))
  }
  ###Cov=0 Solnp###
  tmp <- solnp(p0, fun = objectFunc, eqfun = cov0ConstrFunc, eqB = c(0, 0), LB = paramLowerBound)
  cov0Solnp <- list(tmp)
  #### calculate G2, MaxLogL, AIC####
  getValues <- function(solnpList) {
    for (k in length(solnpList):1) {
      if (length(solnpList) == r - 1) {
        modelName <- paste0("ME", k)
        paramSize <- r^2 - 1 - k
      } else {
        modelName <- "Cov=0"
        paramSize <- r^2 - 2
      }

      G2 <- 2 * ((-fullModel(freq)) - (-solnpList[[k]]$value[length(solnpList[[k]]$value)]))
      maxLogLikeli <- moleculeOfConst - denominatorOfConst + (-solnpList[[k]]$value[length(solnpList[[k]]$value)])
      AIC <- -2 * maxLogLikeli + 2 * paramSize

      fittingValue <- round(sum(freq) * (solnpList[[k]]$pars), 3)
      resultMatrix <- t(matrix(paste0(freq, " (", fittingValue, ")"), r, r))

      analysResults <<- append(analysResults, list(list(deviance = G2, df.residual = k, aic = AIC, result = resultMatrix)))
      names(analysResults)[length(analysResults)] <<- modelName
    }
  }

  getValues(MEkSolnp)
  getValues(cov0Solnp)


  ##### show results #####
  anothNames <- dfs <- G2s <- AICs <- pValues <- codes <- c()

  anothNameTargetModels <- paste0(c("LSI", "LSU", "LSQI", "LSQU", "LS", "ME"), r - 1)
  anothNameModels <- paste0("(", c("I", "U", "QI", "QU", "QS", "MH"), ")")

  for (model in analysResults) {
    modelName <- names(analysResults[match(list(model), analysResults)])
    anothName <- ""
    for (i in 1:length(anothNameTargetModels)) {
      if (modelName == anothNameTargetModels[i]) anothName <- anothNameModels[i]
    }
    pValue <- round(1 - pchisq(model$deviance, model$df.residual), 4)
    code <- ""
    if (0.05 < pValue && pValue < 0.1) {
      code <- "."
    } else {
      for (alpha in c(0.05, 0.01, 0.001)) {
        if (pValue < alpha) code <- paste0(code, "*")
      }
    }

    anothNames <- append(anothNames, anothName)
    dfs <- append(dfs, model$df.residual)
    G2s <- append(G2s, round(model$deviance, 2))
    AICs <- append(AICs, round(model$aic, 2))
    pValues <- append(pValues, pValue)
    codes <- append(codes, code)
  }
  resultForDisplay <- data.frame(model = names(analysResults), anothName = anothNames, df = dfs, G2 = G2s, AIC = AICs, pValue = pValues, code = codes)
  names(resultForDisplay)[2] <- names(resultForDisplay)[7] <- ""
  names(resultForDisplay)[4] <- "G^2"
  names(resultForDisplay)[6] <- "Pr(>G^2)"

  globalAnalysResults <<- analysResults
  cat("\n")

  if (is.character(sort)) {
    print(resultForDisplay[order(resultForDisplay[[sort]]), ])
  } else {
    print(resultForDisplay)
  }

  cat("---\n")
  cat("Signif. codes:  0 ?e***?f 0.001 ?e**?f 0.01 ?e*?f 0.05 ?e.?f 0.1 ?e ?f 1\n")
}



detail <- function(model) {
  selectedModelResult <- globalAnalysResults[[model]]
  glmObjSize <- 30
  modelName <- names(globalAnalysResults[match(list(selectedModelResult), globalAnalysResults)])
  cat("Model:", modelName, "\n")

  if (length(selectedModelResult) == glmObjSize) {
    fittingValue <- round(fitted(selectedModelResult), 2)
    resultMatrix <- t(matrix(paste0(inputData, " (", fittingValue, ")"), r, r))

    print(summary(selectedModelResult))
    cat("Data:\n")
    print(resultMatrix)
  } else {
    print(selectedModelResult)
  }
  print("(The parenthesized values are the MLEs of expected frequencies under the selected model)")
}

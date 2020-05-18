
# Mitscherlich yield surface response function

mitschFunc <- function(NtotDose, PtotDose, KtotDose, 
                       Asym, 
                       Rate_N, Rate_P, Rate_K, 
                       Env_N, Env_P, Env_K) {
  Asym * ((1-exp(-Rate_N*(NtotDose+Env_N)))) * ((1-exp(-Rate_P*(PtotDose+Env_P)))) * ((1-exp(-Rate_K*(KtotDose+Env_K))))
} 


# Mitscherlich NPK mixed model

mitsh_nlme <- function(constrain.asym = NULL, constrain.env = NULL, constrain.rate = NULL,
                       startVal, startVector, rhs, data,
                       rateExp = 'exp', # should rate be exp-transformed or stay ordinary?
                       scaleDose = 'ordinary', # should dose be log-transformed or stay ordinary?
                       maxIter = 50, minScale = 1e-8,
                       msVerbose = FALSE) {
  
  if (scaleDose == "log10") {
    data$NtotDose[data$NtotDose == 0] <- 0.001
    data$PtotDose[data$PtotDose == 0] <- 0.001
    data$KtotDose[data$KtotDose == 0] <- 0.001
    data$NtotDose <- log10(data$NtotDose)
    data$PtotDose <- log10(data$PtotDose)
    data$KtotDose <- log10(data$KtotDose)
  }
  
  # Build formula according to specified constrains
  if (is.null(constrain.asym)) {
    asymTag <- 'Asym'
  } else {
    asymTag <- paste(constrain.asym[1], '+', constrain.asym[2], '/(1+exp(-Asym))')
  }
  
  if (is.null(constrain.env)) {
    envTagN = 'Env_N'
    envTagP = 'Env_P'
    envTagK = 'Env_K'
  } else {
    envTagN <- paste(constrain.env[1], '+', constrain.env[2], '/(1+exp(-Env_N))')
    envTagP <- paste(constrain.env[1], '+', constrain.env[2], '/(1+exp(-Env_P))')
    envTagK <- paste(constrain.env[1], '+', constrain.env[2], '/(1+exp(-Env_K))')
  }
  
  if (is.null(constrain.rate)) {
    rateTagN = 'Rate_N' # default is exponential
    rateTagP = 'Rate_P'
    rateTagK = 'Rate_K'
  } else {
    rateTagN <- paste(constrain.rate[1], '+', constrain.rate[2], '/(1+exp(-Rate_N))')
    rateTagP <- paste(constrain.rate[1], '+', constrain.rate[2], '/(1+exp(-Rate_P))')
    rateTagK <- paste(constrain.rate[1], '+', constrain.rate[2], '/(1+exp(-Rate_K))')
  }
  
  if (rateExp == 'exp') {
    mainMod <- as.formula(paste('RendVendable ~', 
                                asymTag, '* ((1 - exp( -exp(', rateTagN, ') * (NtotDose +', envTagN, '))))',
                                         '* ((1 - exp( -exp(', rateTagP, ') * (PtotDose +', envTagP, '))))', 
                                         '* ((1 - exp( -exp(', rateTagK, ') * (KtotDose +', envTagK, '))))'))
  } else {
    mainMod <- as.formula(paste('RendVendable ~', 
                                asymTag, '* ((1 - exp( -', rateTagN, '* (NtotDose +', envTagN, '))))',
                                         '* ((1 - exp( -', rateTagP, '* (NtotDose +', envTagP, '))))',
                                         '* ((1 - exp( -', rateTagK, '* (NtotDose +', envTagK, '))))'))
  }
  
  mm_NPK <- nlme(mainMod,
                 data = data, 
                 start = c(Asym = startVal[1], start_vector,
                           Rate_N = startVal[2], start_vector,
                           Env_N = startVal[3], start_vector,
                           Rate_P = startVal[4], start_vector,
                           Env_P = startVal[5], start_vector,
                           Rate_K = startVal[6], start_vector,
                           Env_K = startVal[7], start_vector
                           ), 
                 fixed = list(as.formula(paste("Asym ~ ", rhs)),
                              as.formula(paste("Rate_N ~ ", rhs)),
                              as.formula(paste("Env_N ~ ", rhs)),
                              as.formula(paste("Rate_P ~ ", rhs)),
                              as.formula(paste("Env_P ~ ", rhs)),
                              as.formula(paste("Rate_K ~ ", rhs)),
                              as.formula(paste("Env_K ~ ", rhs))
                              ), 
                 random = Asym ~ 1 | NoEssai/NoBloc,
                 control = list(maxIter=100, returnObject=TRUE, 
                                msVerbose=FALSE, minScale=1e-8),
                 method = 'REML')
  
  return(mm_NPK)
}


# Yield response prediction from new data, by hand
## Mitscherlich prediction

pred_mitsch <- function(
  mm, 
  newdata, 
  rhs, 
  NtotDose, PtotDose, KtotDose, 
  constrain.asym = NULL, 
  constrain.env = NULL, 
  constrain.rate = NULL, 
  ranEf = 0#, 
  #rateExp = 'exp'
  ) {
  library(stringi)
  # parameters
  ## collect parameters from the model
  parameter <- summary(mm)$tTable
  
  ## collect specific mean and standard error parameters for Asymptote, Taux and Environnement
  ### Asym  
  asymParam <- parameter[str_detect(rownames(parameter), "Asym"), 1:2]
  ### Rate_
  tauxParamN <- parameter[str_detect(rownames(parameter), "Rate_N"), 1:2] 
  tauxParamP <- parameter[str_detect(rownames(parameter), "Rate_P"), 1:2] 
  tauxParamK <- parameter[str_detect(rownames(parameter), "Rate_K"), 1:2] 
  ### Env_
  enviParamN <- parameter[str_detect(rownames(parameter), "Env_N"), 1:2]
  enviParamP <- parameter[str_detect(rownames(parameter), "Env_P"), 1:2]
  enviParamK <- parameter[str_detect(rownames(parameter), "Env_K"), 1:2]
  
  ## the linear combination of site specific condition (model.matrix) used to compute Mitscherlich parameters
  modmat <- model.matrix(as.formula(paste('~', paste(rhs, collapse = '+'))), data = newdata)
  ### Asym
  asymMitsch <- modmat %*% asymParam[, 1]
  ### Rate_
  tauxMitschN <- modmat %*% tauxParamN[, 1]
  tauxMitschP <- modmat %*% tauxParamP[, 1]
  tauxMitschK <- modmat %*% tauxParamK[, 1]
  ### Env_
  enviMitschN <- modmat %*% enviParamN[, 1]
  enviMitschP <- modmat %*% enviParamP[, 1]
  enviMitschK <- modmat %*% enviParamK[, 1]
  
  ## Treatments columns
  NtotDose <- newdata$NtotDose
  PtotDose <- newdata$PtotDose
  KtotDose <- newdata$KtotDose
  # predict
  #doseMitsch <- newdata[col_dose]
  
  # Build formula according to specified constrains
  if (is.null(constrain.asym)) {
    asymTag <- asymMitsch[1] + ranEf
  } else {
    asymTag <- constrain.asym[1] + constrain.asym[2]/(1+exp(-(asymMitsch[1] + ranEf)))
  }
  
  if (is.null(constrain.env)) {
    envTagN <- enviMitschN[1]
    envTagP <- enviMitschP[1]
    envTagK <- enviMitschK[1]
  } else {
    envTagN <- constrain.env[1] + constrain.env[2]/(1+exp(-(enviMitschN[1])))
    envTagP <- constrain.env[1] + constrain.env[2]/(1+exp(-(enviMitschP[1])))
    envTagK <- constrain.env[1] + constrain.env[2]/(1+exp(-(enviMitschK[1])))
  }
  
  if (is.null(constrain.rate)) {
    rateTagN <- tauxMitschN[1]
    rateTagP <- tauxMitschP[1]
    rateTagK <- tauxMitschK[1]
  } else {
    rateTagN <- constrain.rate[1] + constrain.rate[2]/(1+exp(-(tauxMitschN[1])))
    rateTagP <- constrain.rate[1] + constrain.rate[2]/(1+exp(-(tauxMitschP[1])))
    rateTagK <- constrain.rate[1] + constrain.rate[2]/(1+exp(-(tauxMitschK[1])))
  }
  
  pred <- mitschFunc(NtotDose, PtotDose, KtotDose, 
                     asymTag, 
                     rateTagN, rateTagP, rateTagK,
                     envTagN, envTagP, envTagK#, 
                     #rateExp = rateExp
                     )
  
  # output
  #  return(pred = pred[, 1])
  return(list(pred = pred, 
              metaParam = list(Asym = asymMitsch[1], 
                               Rate_N = tauxMitschN[1], 
                               Rate_P = tauxMitschP[1],
                               Rate_K = tauxMitschK[1],
                               Env_N = enviMitschN[1],
                               Env_P = enviMitschP[1],
                               Env_K = enviMitschK[1]
                               ),
              drag = list(Asym = t(t(modmat[1, ] * asymParam[, 1])),
                          Rate_N = t(t(modmat[1, ] * tauxParamN[, 1])),
                          Rate_P = t(t(modmat[1, ] * tauxParamP[, 1])),
                          Rate_K = t(t(modmat[1, ] * tauxParamK[, 1])),
                          Env_N = t(t(modmat[1, ] * enviParamN[, 1])), 
                          Env_P = t(t(modmat[1, ] * enviParamP[, 1])),
                          Env_K = t(t(modmat[1, ] * enviParamK[, 1]))))
         )
}

# Regression analysis
## R-squared

rsq <- function(y, y_hat) {
  sum((y_hat - mean(y))^2) / sum((y - mean(y))^2)
}
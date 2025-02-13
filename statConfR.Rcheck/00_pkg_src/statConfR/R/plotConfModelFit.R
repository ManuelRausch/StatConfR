#' @title Plot the prediction of fitted parameters of one model of confidence over the corresponding  data
#'
#' @description The `plotConfModelFit` function plots the predicted distribution of discrimination responses
#' and confidence ratings created from a `data.frame` of parameters obtaind from \code{\link{fitConfModels}}
#' and overlays the predicted distribution over the data to which the model parameters were fitted.
#'
#' @param data  a `data.frame` where each row is one trial, containing following
#' variables:
#' * \code{diffCond} (optional; different levels of discriminability,
#'    should be a factor with levels ordered from hardest to easiest),
#' * \code{rating} (discrete confidence judgments, should be a factor with levels
#'    ordered from lowest confidence to highest confidence;
#'    otherwise will be transformed to factor with a warning),
#' * \code{stimulus} (stimulus category in a binary choice task,
#'    should be a factor with two levels, otherwise it will be transformed to
#'    a factor with a warning),
#' * \code{correct} (encoding whether the response was correct; should  be 0 for
#'    incorrect responses and 1 for correct responses)
#' * \code{participant} (some group ID, most often a participant identifier;
#'    the models given in the second argument are fitted to each subset of `data`
#'    determined by the different values of this column)
#'
#' @param fitted_pars a `data.frame` with one row for each participant and model parameters in different columns.
#' fitted_pars also may contain a column called `model` specifying the model to be visualized.
#' If there is no model column in data or if there are multiple models in fitted_pars,
#' it is necessary to specify the model argument.
#'
#' @param model `character`. See \code{\link{fitConfModels}} for all available models
#'
#' @return a `ggplot` object with empirically observed distribution of responses and confidence ratings
#'  as bars on the x-axis as a function of discriminability (in the rows) and stimulus
#'  (in the columns). Superimposed on the empirical data,
#'  the plot also shows the prediction of one selected model as dots.
#'
#' @examples
#' # 1. Fit some models to each subject of the masked orientation discrimination experiment
#'   # Normally, the fits should be created using the function fitConfModels
#'   # Fits <- fitConfModels(data, models = "WEV", .parallel = TRUE)
#'   # Here, we create the dataframe manually because fitting models takes about
#'   # 10 minutes per model fit per participant on a 2.8GHz processor.
#'   pars <- data.frame(participant = 1:16,
#'   d_1 = c(0.20, 0.05, 0.41, 0.03, 0.00, 0.01, 0.11, 0.03, 0.19, 0.08, 0.00,
#'   0.24, 0.00, 0.00, 0.25, 0.01),
#'   d_2 = c(0.61, 0.19, 0.86, 0.18, 0.17, 0.39, 0.69, 0.14, 0.45, 0.30, 0.00,
#'   0.27, 0.00, 0.05, 0.57, 0.23),
#'   d_3 = c(1.08, 1.04, 2.71, 2.27, 1.50, 1.21, 1.83, 0.80, 1.06, 0.68, 0.29,
#'   0.83, 0.77, 2.19, 1.93, 0.54),
#'   d_4 = c(3.47, 4.14, 6.92, 4.79, 3.72, 3.24, 4.55, 2.51, 3.78, 2.40, 1.95,
#'   2.55, 4.59, 4.27, 4.08, 1.80),
#'   d_5 = c(4.08, 5.29, 7.99, 5.31, 4.53, 4.66, 6.21, 4.67, 5.85, 3.39, 3.39,
#'   4.42, 6.48, 5.35, 5.28, 2.87),
#'   c = c(-0.30, -0.15, -1.37, 0.17, -0.12, -0.19, -0.12, 0.41, -0.27, 0.00,
#'   -0.19, -0.21, -0.91, -0.26, -0.20, 0.10),
#'   theta_minus.4 = c(-2.07, -2.04, -2.76, -2.32, -2.21, -2.33, -2.27, -2.29,
#'   -2.69, -3.80, -2.83, -1.74, -2.58, -3.09, -2.20, -1.57),
#'   theta_minus.3 = c(-1.25, -1.95, -1.92, -2.07, -1.62, -1.68, -2.04, -2.02,
#'   -1.84, -3.37, -1.89, -1.44, -2.31, -2.08, -1.53, -1.46),
#'   theta_minus.2 = c(-0.42, -1.40, -0.37, -1.96, -1.45, -1.27, -1.98, -1.66,
#'   -1.11, -2.69, -1.60, -1.25, -2.21, -1.68, -1.08, -1.17),
#'   theta_minus.1 = c(0.13, -0.90,  0.93, -1.71, -1.25, -0.59, -1.40, -1.00,
#'   -0.34, -1.65, -1.21, -0.76, -1.99, -0.92, -0.28, -0.99),
#'   theta_plus.1 = c(-0.62, 0.82, -2.77, 2.01, 1.39, 0.60, 1.51, 0.90, 0.18,
#'   1.62, 0.99,0.88, 1.67, 0.92, 0.18,  0.88),
#'   theta_plus.2 = c(0.15, 1.45, -1.13,2.17, 1.61, 1.24, 1.99, 1.55, 0.96, 2.44,
#'   1.53, 1.66, 2.00, 1.51, 1.08, 1.05),
#'   theta_plus.3 = c(1.40, 2.24, 0.77, 2.32, 1.80, 1.58, 2.19, 2.19, 1.54, 3.17,
#'   1.86, 1.85, 2.16, 2.09, 1.47, 1.70),
#'   theta_plus.4 = c(2.19, 2.40, 1.75, 2.58, 2.53, 2.24, 2.59, 2.55, 2.58, 3.85,
#'   2.87, 2.15, 2.51, 3.31, 2.27, 1.79),
#'   sigma = c(1.01, 0.64, 1.33, 0.39, 0.30, 0.75, 0.75, 1.07, 0.65, 0.29, 0.31,
#'   0.78, 0.39, 0.42, 0.69, 0.52),
#'   w = c(0.54, 0.50, 0.38, 0.38, 0.36, 0.44, 0.48, 0.48, 0.52, 0.46, 0.53, 0.48,
#'    0.29, 0.45, 0.51, 0.63))
#'
#' # 2. Plot the predicted probabilities based on model and fitted parameters
#'   # against the observed relative frequencies.
#'
#'   PlotFitWEV <- plotConfModelFit(MaskOri, pars, model="WEV")
#'   PlotFitWEV
#'
#' @import ggplot2
#' @importFrom plyr ddply summarise
#' @importFrom Rmisc summarySEwithin
#' @importFrom stats aggregate

#' @author Manuel Rausch, \email{manuel.rausch@hochschule-rhein-waal.de}

#' @export
plotConfModelFit <- function(data, fitted_pars, model = NULL){
  if(is.null(model)){
    if("model" %in% colnames(fitted_pars)){
      if(length(unique(fitted_pars$model))==1){
        model <- unique(fitted_pars$model)
      } else {
        stop("Please use the model argument to specify which model should be used")
      }
    }else {
      stop("Please specify which model should be used")
    }
  }
  if (is.null(data$diffCond)) data$diffCond <- factor(1)
  if (!is.factor(data$diffCond)) {
    data$diffCond <- factor(data$diffCond)
  }
  if(length(unique(data$stimulus)) != 2) {
    stop("There must be exacltly two different possible values of stimulus")
  }

  if (!is.factor(data$stimulus)) {
    data$stimulus <- factor(data$stimulus)
  }
  if (!is.factor(data$rating)) {
    data$rating <- factor(data$rating)
  }
  if(!all(data$correct %in% c(0,1))) stop("correct should be 1 or 0")

  PlotName <-
    switch(model,
           'GN' = "Gaussian noise model",
           'IG' = "Independent Gaussian model",
           'ITGc'  = "Independent truncated Gaussian model: HMetad-Version",
           'ITGcm' = "Independent truncated Gaussian model: Meta-d'-Version",
           'logN' = "Logistic noise model",
           'logWEV' = "Logistic weighted evidence and visibility model",
           'PDA' = "Post-decisional accumulation model",
           'WEV' = "Weighted evidence and visibility model",
           'SDT' = "Signal detection rating model") # models are color coded

  # 1. First aggregate on the level of subjects

  AggDist <-
    plyr::ddply(data,
                ~  diffCond * rating *
                  stimulus * correct * participant, #,
                plyr::summarise,
                p = length(.data$rating),  .drop=FALSE)

  AggDist <- plyr::ddply(AggDist, ~
                           diffCond * stimulus * participant,
                         transform, N = sum(.data$p))
  AggDist$p <- AggDist$p / AggDist$N


  AggDist$rating <- as.numeric(AggDist$rating)
  AggDist$correct <-
    factor(AggDist$correct)
  levels(AggDist$correct)  <-
    c("incorrect", "correct")

  # 2. aggregate across subjects
  AggDist <-
    Rmisc::summarySEwithin(AggDist , measurevar = "p",
                           withinvars = c("diffCond", "correct", "rating", "stimulus"),
                           idvar = "participant",
                           na.rm = TRUE, .drop = TRUE)
  AggDist$rating <- as.numeric(AggDist$rating)
  levels(AggDist$stimulus) <- c("S = -1", "S = 1")
  AggDist$diffCond <-
    factor(as.numeric(AggDist$diffCond)) # diffCond should code the order of difficulty levels
  levels(AggDist$diffCond) <-
    paste("K =", as.numeric(levels(AggDist$diffCond)))

  # 4) create the prediction from model fit

  predictDataFun <- switch(model,
                           'WEV' = predictDataWEV,
                           'logWEV' = predictDataLogWEV,
                           'SDT' = predictDataSDT,
                           'IG' = predictData2Chan,
                           'ITGc' = predictDataIndTruncF,
                           'ITGcm' = predictDataIndTruncML,
                           'GN' = predictDataNoisy,
                           'PDA' = predictDataISDT,
                           'logN' = predictDataLognorm)

  if("model" %in% colnames(fitted_pars)){
    if(length(unique(fitted_pars$model))>1){
      fitted_pars <- fitted_pars[fitted_pars$model == model,]
    }}

  predictedData <- plyr::ddply(fitted_pars, ~participant, predictDataFun)
  predictedData$correct <-
    factor(ifelse(predictedData$stimulus==predictedData$response, "correct", "incorrect"))

  AggPredictedData <-
    aggregate(p ~ stimulus*diffCond*rating*correct, predictedData, mean)
  AggPredictedData$stimulus <- factor(AggPredictedData$stimulus)
  levels(AggPredictedData$stimulus) <- c("S = -1", "S = 1")
  AggPredictedData$diffCond <- factor(AggPredictedData$diffCond)
  levels(AggPredictedData$diffCond) <-
    paste("K =", as.numeric(levels(AggPredictedData$diffCond)))


  # 5)  create a plot with the observed data
  PlotObsVsPred  <-
    ggplot(AggDist,
           aes(x = .data$rating, y = .data$p)) +
    facet_grid(diffCond ~ stimulus + correct) +  # No .data$ needed
    geom_bar(stat = "identity",
             fill = "white", color = "black") +
    geom_errorbar(aes(ymin = .data$p - .data$se, ymax = .data$p + .data$se), width = 0) +
    geom_point(data = AggPredictedData, aes(x = .data$rating, y = .data$p)) + # Fix
    xlab("Confidence rating") +
    ylab("Probability") +
    coord_cartesian(ylim = c(0,1)) +  # Avoids ggplot removing data
    ggtitle(PlotName) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme_minimal()

  # create a plot from the observed data
  PlotObsVsPred
}

# internal functions to generate the predicted probabilites
predictDataSDT <- function(paramDf){

  nCond <- sum(is.finite(c(t(paramDf[,grep(pattern = "d_", names(paramDf), value=T)]))))
  nRatings <- sum(is.finite(c(t(paramDf[,grep(pattern = "^theta_minus.", names(paramDf), value=T)]))))+1
  ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
  locA <- -ds/2
  locB <- ds/2
  theta <- paramDf$c
  c_RA <- c(-Inf,c(t(paramDf[,paste("theta_minus.", (nRatings-1):1, sep="")])), theta)
  c_RB <- c(theta,c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])), Inf)

  p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

  P_SBRB <- Vectorize(function(j,i) pnorm(q=c_RB[i+1], locB[j]) - pnorm(q=c_RB[i], locB[j]))
  P_SBRA <- Vectorize(function(j,i) pnorm(q=c_RA[i+1], locB[j]) - pnorm(q=c_RA[i], locB[j]))
  P_SARA <- Vectorize(function(j,i) pnorm(q=c_RA[i+1], locA[j]) - pnorm(q=c_RA[i], locA[j]))
  P_SARB <- Vectorize(function(j,i) pnorm(q=c_RB[i+1], locA[j]) - pnorm(q=c_RB[i], locA[j]))

  p_SB_RB <- cbind(stimulus = 1, response = 1,
                   plyr::mdply(.data=p_SB_RB, .fun= P_SBRB))
  p_SB_RA <- cbind(stimulus = 1, response = -1,
                   plyr::mdply(.data=p_SB_RA, .fun= P_SBRA))
  p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
  p_SA_RA <- cbind(stimulus = -1, response = -1,
                   plyr::mdply(.data=p_SA_RA, .fun= P_SARA))
  p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
  p_SA_RB <- cbind(stimulus = -1, response = 1,
                   plyr::mdply(.data=p_SA_RB, .fun= P_SARB))

  res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
  colnames(res) <- c("stimulus", "response", "diffCond",
                     "rating", "p")
  res$p[is.na(res$p) | is.nan(res$p)] <- 0
  res
}

# (ii) Noisy

predictDataNoisy <-
  function(paramDf){

    nCond <- length(grep(pattern = "d_", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "^theta_minus.", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
    locA <- -ds/2
    locB <- ds/2
    theta <- paramDf$c
    c_RA <- c(-Inf,c(t(paramDf[,paste("theta_minus.", (nRatings-1):1, sep="")])), Inf)
    c_RB <- c(-Inf,c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])), Inf)
    sigma <- paramDf$sigma

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

    P_SBRB_Noisy <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locB[j]) * (pnorm(q=c_RB[i+1], x, sigma) - pnorm(q=c_RB[i], x, sigma)),
                lower = theta,
                upper = Inf,
                rel.tol = 10^-8)$value
    })
    P_SBRA_Noisy <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locB[j]) * (pnorm(q=c_RA[i+1], x, sigma) - pnorm(q=c_RA[i], x, sigma)),
                lower = -Inf,
                upper = theta,
                rel.tol = 10^-8)$value
    })
    P_SARA_Noisy <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locA[j]) * (pnorm(q=c_RA[i+1], x, sigma) - pnorm(q=c_RA[i], x, sigma)),
                lower = -Inf,
                upper = theta,
                rel.tol = 10^-8)$value
    })
    P_SARB_Noisy <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locA[j]) * (pnorm(q=c_RB[i+1], x, sigma) - pnorm(q=c_RB[i], x, sigma)),
                lower = theta,
                upper= Inf,
                rel.tol = 10^-8)$value
    })

    p_SB_RB <- cbind(stimulus = 1, response = 1,
                     plyr::mdply(.data=p_SB_RB, .fun= P_SBRB_Noisy))
    p_SB_RA <- cbind(stimulus = 1, response = -1,
                     plyr::mdply(.data=p_SB_RA, .fun= P_SBRA_Noisy))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = -1, response = -1,
                     plyr::mdply(.data=p_SA_RA, .fun= P_SARA_Noisy))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = -1, response = 1,
                     plyr::mdply(.data=p_SA_RB, .fun= P_SARB_Noisy))

    res <- rbind(p_SB_RB,p_SB_RA,  p_SA_RA,p_SA_RB)
    colnames(res) <- c("stimulus", "response", "diffCond", "rating", "p")
    res$p[is.na(res$p) | is.nan(res$p)] <- 0
    res
  }


# (iii) PDA

predictDataISDT <-
  function(paramDf){

    nCond <- length(grep(pattern = "^d_", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "^theta_minus.", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
    locA <- -ds/2
    locB <- ds/2
    theta <- paramDf$c
    a <-  paramDf$b
    sigma <- sqrt(a)

    c_RA <- c(-Inf, c(t(paramDf[,paste("theta_minus.", (nRatings-1):1, sep="")])), Inf)
    c_RB <- c(-Inf, c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])), Inf)

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

    P_SBRB <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locB[j]) * (pnorm(q=c_RB[i+1], x + locB[j]*a, sigma) - pnorm(q=c_RB[i], x + locB[j]*a, sigma)),
                lower = theta,
                upper = Inf,
                rel.tol = 10^-8)$value
    })
    P_SBRA <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locB[j]) * (pnorm(q=c_RA[i+1], x + locB[j]*a, sigma) - pnorm(q=c_RA[i], x + locB[j]*a, sigma)),
                lower = -Inf,
                upper = theta,
                rel.tol = 10^-8)$value
    })
    P_SARA <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locA[j]) * (pnorm(q=c_RA[i+1], x + locA[j]*a, sigma) - pnorm(q=c_RA[i], x + locA[j]*a, sigma)),
                lower = -Inf,
                upper = theta,
                rel.tol = 10^-8)$value
    })
    P_SARB <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locA[j]) * (pnorm(q=c_RB[i+1], x + locA[j]*a, sigma) - pnorm(q=c_RB[i], x +  locA[j]*a, sigma)),
                lower = theta,
                upper= Inf,
                rel.tol = 10^-8)$value
    })

    p_SB_RB <- cbind(stimulus = 1, response = 1,
                     plyr::mdply(.data=p_SB_RB, .fun= P_SBRB))
    p_SB_RA <- cbind(stimulus = 1, response = -1,
                     plyr::mdply(.data=p_SB_RA, .fun= P_SBRA))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = -1, response = -1,
                     plyr::mdply(.data=p_SA_RA, .fun= P_SARA))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = -1, response = 1,
                     plyr::mdply(.data=p_SA_RB, .fun= P_SARB))

    res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
    colnames(res) <- c("stimulus", "response", "diffCond",
                       "rating", "p")
    res$p[is.na(res$p) | is.nan(res$p)] <- 0
    res
  }

# (iv) IG

predictData2Chan <- function(paramDf){
  nCond <- sum(is.finite(c(t(paramDf[,grep(pattern = "d_", names(paramDf), value=T)]))))
  nRatings <- sum(is.finite(c(t(paramDf[,grep(pattern = "^theta_minus.", names(paramDf), value=T)]))))+1
  ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
  locA1 <- -ds/2
  locB1 <- ds/2
  theta <- paramDf$c
  c_RA <- c(-Inf,c(t(paramDf[,paste("theta_minus.", (nRatings-1):1, sep="")])), Inf)
  c_RB <- c(-Inf,c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])), Inf)
  metads <- paramDf$m * ds
  locA2 <- -metads/2
  locB2 <- metads/2

  p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

  P_SBRB <- Vectorize(function(j,i){
    (1 - pnorm(theta, locB1[j])) * (pnorm(c_RB[i+1], locB2[j]) - pnorm(c_RB[i], locB2[j]))
  })
  P_SBRA <- Vectorize(function(j,i){
    pnorm(theta, locB1[j]) * (pnorm(c_RA[i+1], locB2[j]) - pnorm(c_RA[i], locB2[j]))
  })
  P_SARA <-  Vectorize(function(j,i){
    pnorm(theta, locA1[j]) * (pnorm(c_RA[i+1], locA2[j]) - pnorm(c_RA[i], locA2[j]))
  })
  P_SARB <-  Vectorize(function(j,i){
    (1 - pnorm(theta, locA1[j])) * (pnorm(c_RB[i+1], locA2[j]) - pnorm(c_RB[i], locA2[j]))
  })

  p_SB_RB <- cbind(stimulus = 1, response = 1,
                   plyr::mdply(.data=p_SB_RB, .fun= P_SBRB))
  p_SB_RA <- cbind(stimulus = 1, response = -1,
                   plyr::mdply(.data=p_SB_RA, .fun= P_SBRA))
  p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
  p_SA_RA <- cbind(stimulus = -1, response = -1,
                   plyr::mdply(.data=p_SA_RA, .fun= P_SARA))
  p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
  p_SA_RB <- cbind(stimulus = -1, response = 1,
                   plyr::mdply(.data=p_SA_RB, .fun= P_SARB))

  res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
  colnames(res) <- c("stimulus", "response", "diffCond",
                     "rating", "p")
  res$p[is.na(res$p) | is.nan(res$p)] <- 0
  res
}

# (v) WEV

predictDataWEV <-
  function(paramDf){

    nCond <- length(grep(pattern = "d_", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "^theta_minus.", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
    locA <- -ds/2
    locB <- ds/2
    c_RA <- c(-Inf,c(t(paramDf[,paste("theta_minus.", (nRatings-1):1, sep="")])), Inf)
    c_RB <- c(-Inf,c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])), Inf)

    theta <- paramDf$c
    w <- paramDf$w
    sigma <- paramDf$sigma

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

    P_SBRB_CEV <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j]) * (pnorm(q=c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - pnorm(q=c_RB[i], (1 - w) *x + ds[j] * w, sigma)),
        lower = theta,
        upper = Inf,
        rel.tol = 10^-8)$value
    })
    P_SBRA_CEV <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j]) * (pnorm(q=c_RA[i+1], (1 - w) * x - ds[j] * w, sigma) - pnorm(q=c_RA[i], (1 - w) *x - ds[j] * w, sigma)),
        lower = -Inf,
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARA_CEV <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j]) * (pnorm(q=c_RA[i+1], (1 - w) * x - ds[j] * w, sigma) - pnorm(q=c_RA[i], (1 - w) *x - ds[j] * w, sigma)),
        lower = -Inf,
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARB_CEV <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j]) * (pnorm(q=c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - pnorm(q=c_RB[i], (1 - w) *x + ds[j] * w, sigma)),
        lower = theta,
        upper= Inf,
        rel.tol = 10^-8)$value
    })

    p_SB_RB <- cbind(stimulus = 1, response = 1,
                     plyr::mdply(.data=p_SB_RB, .fun= P_SBRB_CEV))
    p_SB_RA <- cbind(stimulus = 1, response = -1,
                     plyr::mdply(.data=p_SB_RA, .fun= P_SBRA_CEV))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = -1, response = -1,
                     plyr::mdply(.data=p_SA_RA, .fun= P_SARA_CEV))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = -1, response = 1,
                     plyr::mdply(.data=p_SA_RB, .fun= P_SARB_CEV))

    res <- rbind(p_SB_RB,p_SB_RA,  p_SA_RA,p_SA_RB)
    colnames(res) <- c("stimulus", "response", "diffCond", "rating", "p")
    res$p[is.na(res$p) | is.nan(res$p)] <- 0
    res
  }


# (vi) ITGcm

predictDataIndTruncML <-
  function(paramDf){

    nCond <- length(grep(pattern = "d_", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "^theta_minus.", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
    locA1 <- -ds/2
    locB1 <- ds/2
    theta <- paramDf$c

    m_ratio <- paramDf$m
    metads <- m_ratio * ds
    locA2 <- -metads/2
    locB2 <- metads/2
    meta_c <- m_ratio * theta

    c_RA <- c(-Inf, c(t(paramDf[,paste("theta_minus.", (nRatings-1):1, sep="")])), meta_c)
    c_RB <- c(meta_c, c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])), Inf)

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

    P_SBRB <- Vectorize(function(j,i){
      (1 - pnorm(theta, locB1[j])) *
        (pnorm(c_RB[i+1], locB2[j]) - pnorm(c_RB[i], locB2[j])) /
        (1 - pnorm(meta_c, locB2[j]))
    })
    P_SBRA <- Vectorize(function(j,i){
      pnorm(theta, locB1[j]) *
        (pnorm(c_RA[i+1], locB2[j]) - pnorm(c_RA[i], locB2[j])) /
        pnorm(meta_c, locB2[j])
    })
    P_SARA <-  Vectorize(function(j,i){
      pnorm(theta, locA1[j]) *
        (pnorm(c_RA[i+1], locA2[j]) - pnorm(c_RA[i], locA2[j])) /
        pnorm(meta_c, locA2[j])
    })
    P_SARB <-  Vectorize(function(j,i){
      (1 - pnorm(theta, locA1[j])) *
        (pnorm(c_RB[i+1], locA2[j]) - pnorm(c_RB[i], locA2[j])) /
        (1 - pnorm(meta_c, locA2[j]))
    })

    p_SB_RB <- cbind(stimulus = 1, response = 1,
                     plyr::mdply(.data=p_SB_RB, .fun= P_SBRB))
    p_SB_RA <- cbind(stimulus = 1, response = -1,
                     plyr::mdply(.data=p_SB_RA, .fun= P_SBRA))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = -1, response = -1,
                     plyr::mdply(.data=p_SA_RA, .fun= P_SARA))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = -1, response = 1,
                     plyr::mdply(.data=p_SA_RB, .fun= P_SARB))

    res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
    colnames(res) <- c("stimulus", "response", "diffCond",
                       "rating", "p")
    res$p[is.na(res$p)] <- 0
    res$p[is.nan(res$p)] <- 0
    res$p[res$p < 0] <- 0

    res
  }


# (vii) ITGc

predictDataIndTruncF <-
  function(paramDf){

    nCond <- length(grep(pattern = "d_", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "^theta_minus.", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
    locA1 <- -ds/2
    locB1 <- ds/2
    theta <- paramDf$c

    m_ratio <- paramDf$m
    metads <- m_ratio * ds
    locA2 <- -metads/2
    locB2 <- metads/2
    meta_c <-  theta

    c_RA <- c(-Inf, c(t(paramDf[,paste("theta_minus.", (nRatings-1):1, sep="")])), meta_c)
    c_RB <- c(meta_c, c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])), Inf)

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

    P_SBRB <- Vectorize(function(j,i){
      (1 - pnorm(theta, locB1[j])) * (pnorm(c_RB[i+1], locB2[j]) - pnorm(c_RB[i], locB2[j])) / (1 - pnorm(meta_c, locB2[j]))
    })
    P_SBRA <- Vectorize(function(j,i){
      pnorm(theta, locB1[j]) * (pnorm(c_RA[i+1], locB2[j]) - pnorm(c_RA[i], locB2[j])) / pnorm(meta_c, locB2[j])
    })
    P_SARA <-  Vectorize(function(j,i){
      pnorm(theta, locA1[j]) *
        (pnorm(c_RA[i+1], locA2[j]) - pnorm(c_RA[i], locA2[j])) /
        pnorm(meta_c, locA2[j])
    })
    P_SARB <-  Vectorize(function(j,i){
      (1 - pnorm(theta, locA1[j])) *
        (pnorm(c_RB[i+1], locA2[j]) - pnorm(c_RB[i], locA2[j])) /
        (1 - pnorm(meta_c, locA2[j]))
    })

    p_SB_RB <- cbind(stimulus = 1, response = 1,
                     plyr::mdply(.data=p_SB_RB, .fun= P_SBRB))
    p_SB_RA <- cbind(stimulus = 1, response = -1,
                     plyr::mdply(.data=p_SB_RA, .fun= P_SBRA))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = -1, response = -1,
                     plyr::mdply(.data=p_SA_RA, .fun= P_SARA))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = -1, response = 1,
                     plyr::mdply(.data=p_SA_RB, .fun= P_SARB))

    res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
    colnames(res) <- c("stimulus", "response", "diffCond",
                       "rating", "p")
    res$p[is.na(res$p)] <- 0
    res$p[is.nan(res$p)] <- 0
    res$p[res$p < 0] <- 0
    res
  }

# (viii) logN

predictDataLognorm <- function(paramDf){
  nCond <- length(grep(pattern = "d_", names(paramDf), value=T))
  nRatings <- (length(grep(pattern = "^M_theta_minus.", names(paramDf), value=T)))+1
  ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
  locA <- -ds/2
  locB <- ds/2
  theta <- paramDf$c
  sigma <- paramDf$sigma

  loc_RA <-  c(Inf,log(abs(c(t(paramDf[,paste("M_theta_minus.", (nRatings-1):1, sep="")])) - theta)) -
                 .5*sigma^2, -Inf)  # the order here is REVERSED!!!
  loc_RB <- c(-Inf, log(c(t(paramDf[,paste("M_theta_plus.", 1:(nRatings-1), sep="")])) - theta) -
                .5*sigma^2, Inf)

  p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

  P_SBRB <-  Vectorize(function(j,i){
    integrate(function(x) dnorm(x, locB[j]) *
                (plnorm(x-theta, loc_RB[i], sigma) - plnorm(x-theta, loc_RB[i+1], sigma)),
              lower = theta, upper = Inf, rel.tol = 10^-8)$value
  })
  P_SBRA <-  Vectorize(function(j,i){
    integrate(function(x) dnorm(x, locB[j]) *
                (plnorm(theta-x,  loc_RA[i+1], sigma) - plnorm(theta-x,  loc_RA[i], sigma)),
              lower = -Inf, upper =theta, rel.tol = 10^-8)$value
  })
  P_SARA <- Vectorize(function(j,i){
    integrate(function(x) dnorm(x, locA[j]) *
                (plnorm(theta-x,  loc_RA[i+1], sigma) - plnorm(theta-x,  loc_RA[i], sigma)),
              lower = -Inf, upper =theta, rel.tol = 10^-8)$value
  })
  P_SARB <- Vectorize(function(j,i){
    integrate(function(x) dnorm(x, locA[j]) *
                (plnorm(x-theta, loc_RB[i], sigma) - plnorm(x-theta, loc_RB[i+1], sigma)),
              lower = theta, upper = Inf, rel.tol = 10^-8)$value
  })

  p_SB_RB <- cbind(stimulus = 1, response = 1,
                   plyr::mdply(.data= expand.grid(j = 1:nCond, i = 1:nRatings), .fun= P_SBRB))
  p_SB_RA <- cbind(stimulus = 1, response = -1,
                   plyr::mdply(.data=expand.grid(j = 1:nCond, i = 1:nRatings), .fun= P_SBRA))
  p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
  p_SA_RA <- cbind(stimulus = -1, response =-1,
                   plyr::mdply(.data=expand.grid(j = 1:nCond, i = 1:nRatings), .fun= P_SARA))
  p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
  p_SA_RB <- cbind(stimulus = -1, response = 1,
                   plyr::mdply(.data=expand.grid(j = 1:nCond, i = 1:nRatings), .fun= P_SARB))

  res <- rbind(p_SB_RB,p_SB_RA,  p_SA_RA,p_SA_RB)
  colnames(res) <- c("stimulus", "response", "diffCond", "rating", "p")
  res$p[is.na(res$p) | is.nan(res$p)] <- 0
  res

}

# (ix) logWEV

predictDataLogWEV <- function(paramDf){
  nCond <- length(grep(pattern = "d_", names(paramDf), value=T))
  nRatings <- (length(grep(pattern = "^theta.minus", names(paramDf), value=T)))+1
  ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
  locA <- -ds/2
  locB <- ds/2
  theta <- paramDf$c
  sigma <- paramDf$sigma
  w <- paramDf$w
  c_RA <- c(0,c(t(paramDf[,paste("theta_minus.", 1:(nRatings-1), sep="")])),
            -Inf) # due to the lognormal distributions, the rating criteria are bounded by 0
  c_RB <- c(0, c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])),
            Inf)


  p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

  P_SBRB <-  Vectorize(function(j,i){
    integrate(
      function(x) dnorm(x, locB[j]) * (plnorm(q = c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - plnorm(q = c_RB[i], (1 - w) * x + ds[j] * w,  sigma)),
      lower = theta,
      upper = Inf,
      rel.tol = 10^-8)$value
  })
  P_SBRA <-  Vectorize(function(j,i){
    integrate(
      function(x) dnorm(x, locB[j]) * (plnorm(q = -c_RA[i+1], -(1 - w) * x + ds[j] * w, sigma) - plnorm(q = -c_RA[i], -(1 - w) * x + ds[j] * w, sigma)),
      lower = -Inf,
      upper = theta,
      rel.tol = 10^-8)$value
  })
  P_SARA <- Vectorize(function(j,i){
    integrate(
      function(x) dnorm(x, locA[j]) * (plnorm(q= -c_RA[i+1], -(1 - w) * x + ds[j] * w, sigma) - plnorm(q = -c_RA[i], -(1 - w) * x + ds[j] * w, sigma)),
      lower = -Inf,
      upper = theta,
      rel.tol = 10^-8)$value
  })
  P_SARB <- Vectorize(function(j,i){
    integrate(
      function(x) dnorm(x, locA[j]) * (plnorm(q=c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - plnorm(q=c_RB[i], (1 - w) * x + ds[j] * w, sigma)),
      lower = theta,
      upper= Inf,
      rel.tol = 10^-8)$value
  })


  p_SB_RB <- cbind(stimulus = 1, response = 1,
                   plyr::mdply(.data= expand.grid(j = 1:nCond, i = 1:nRatings), .fun= P_SBRB))
  p_SB_RA <- cbind(stimulus = 1, response = -1,
                   plyr::mdply(.data=expand.grid(j = 1:nCond, i = 1:nRatings), .fun= P_SBRA))
  p_SA_RA <- cbind(stimulus = -1, response =-1,
                   plyr::mdply(.data=expand.grid(j = 1:nCond, i = 1:nRatings), .fun= P_SARA))
  p_SA_RB <- cbind(stimulus = -1, response = 1,
                   plyr::mdply(.data=expand.grid(j = 1:nCond, i = 1:nRatings), .fun= P_SARB))

  res <- rbind(p_SB_RB,p_SB_RA,  p_SA_RA,p_SA_RB)
  colnames(res) <- c("stimulus", "response", "diffCond", "rating", "p")
  res$p[is.na(res$p) | is.nan(res$p)] <- 0
  res
}

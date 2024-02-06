#' Simulate data  according to a static model of confidence
#'
#' @param model `character` of length 1.
#' Models implemented so far: 'WEV', 'SDT', 'GN', 'PDA', 'IG', 'ITGc', 'ITGcm', 'logN', and 'logWEV'.
#' @param data  a `data.frame` that contains all parameters to simulate a data set,
#' with one row and the different parameters in different columns. Which paramters are needed depends on the specific model:
#' * \code{N} the number of trials be simulated
#' * \code{d_1}, \code{d_2}, ... sensitivity parameters. The number of sensitivity parameters determines the number of conditions
#' * \code{c}: discrimination bias
#' * \code{theta_minus.1}, \code{theta_minus.2},... confidence criteria associated with the response R = -1. The function simulates one more confidence category than there are confidence criteria.
#' * \code{theta_plus.1}, \code{theta_plus.2},... confidence criteria associated with the response R = 1. The function simulates one more confidence category than there are confidence criteria.
#' * \code{w} only for WEV: the visibility weighting parameter, bounded between 0 and 1
#' * \code{sigma} only for WEV, GN, logN, and logWEV: confidence noise, bounded between 0 and Inf.
#' * \code{m} only for IG, ITGm, and ITGcm: metacognitive efficiency parameter, bounded between 0 and Inf
#' * \code{b} only for PDA: postdecisional accumulation parameter, bounded between 0 and Inf
#' * \code{M_theta_minus.1}, \code{M_theta_minus.2},.... Only for logN: Mean confidence criteria associated with the response R = -1.
#' * \code{M_theta_plus.1}, \code{M_theta_plus.2},... . Only for logN: Mean confidence criteria associated with the response R = 1.
#'
#' @return a dataframe with \code{N} rows, and the columns
#' \code{stimulus}, \code{correct} and \code{rating}.
#'
#' @details see \code{fitConf} for a detailed description of the different models.
#'
#' @md
#'
#'
#' @author Manuel Rausch, \email{manuel.rausch@hochschule-rhein-waal.de}

#' @name simConf


#' @importFrom plyr mdply ddply
#'
#' @examples
#' 1. define some parameters
#' paramDf <- data.frame(d_1 = 0, d_2 = 2, d_3 = 4,c = .0,  theta_minus.2 = -2, theta_minus.1 = -1, theta_plus.1 = 1, theta_plus.2 = 2,  sigma = 1/2, w = 0.5, nTrials = 500)
#' 2. Simulate dataset
#' SimulatedData <- simConf(model = "WEV", paramDf)
#'
#' @export
simConf <- function(model = "SDT",  paramDf) {
  AllModels <-
    c('WEV', 'SDT','IG','ITGc', 'ITGcm', 'GN', 'PDA', 'logN', 'logWEV')
  if (! model %in% AllModels){
    stop(paste(model, " not implemented!"))
  }
  if(!"N" %in% colnames(paramDf)){
   stop("Please specify the number of trials in the column N")
  }

  SimFun <- switch(model,
                   'WEV' = generateDataWEV,
                   'logWEV' = generateDataLogWEV,
                   'SDT' = generateDataSDT,
                   'IG' = generateData2Chan,
                   'ITGc' = generateDataIndTruncF,
                   'ITGcm' = generateDataIndTruncML,
                   'GN' = generateDataNoisy,
                   'PDA' = generateDataISDT,
                   'logN' = generateDataLognorm)
  paramDf$V1 <- 1:nrow(paramDf)
  SimData <- plyr::ddply(paramDf, ~V1, SimFun)
  SimData
}




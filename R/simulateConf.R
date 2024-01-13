#' Simulate data  according to a static model of confidence
#'
#' @param model `character` of length 1.
#' Models implemented so far: 'WEV', 'SDT', 'GN', 'PDA', 'IG', 'ITGc', 'ITGcm', 'logN', and 'logWEV'.
#' @param data  a `data.frame` that contains all parameters to simulate a data set,
#' with one row and the different parameters in different columns:
#' * \code{nTirals} the number of trials be simulated
#' * \code{d_1}, \code{d_2}, ... sensitivity parameters, one for each condition.
#' * \code{theta_minus.1}, \code{theta_minus.2},... confidence criteria associated with the response R = -1. The function simulates one more confidence category than there are confidence criteria.
#' * \code{theta_plus.1}, \code{theta_plus.2},... confidence criteria associated with the response R = 1. The function simulates one more confidence category than there are confidence criteria.
#' * \code{w} only for WEV: the visibility weighting parameter, bounded between 0 and 1
#' * \code{sigma} only for WEV, GN, logN, and logWEV: confidecne noise, bounded between 0 and Inf.

#' @return
#'
#' @details seee \code{fitConf} for a detailed description of the different models.

#' @importFrom plyr mdply
#' @export
simConf <- function(model = "SDT",  paramDf) {
  AllModels <-
    c('WEV', 'SDT','IG','ITGc', 'ITGcm', 'GN', 'PDA', 'logN')
  if (! model %in% AllModels){
    stop(paste(paste(setdiff(models, AllModels),collapse = " and "), " not implemented!"))
  }
  SimFun <- switch(model,
                   'WEV' = generateDataWEV,
                   'SDT' = generateDataSDT,
                   'IG' = generateData2Chan,
                   'ITGc' = generateDataIndTruncF,
                   'ITGcm' = generateDataIndTruncML,
                   'GN' = generateDataNoisy,
                   'PDA' = generateDataISDT,
                   'logN' = generateDataLognorm)
  SimData <- SimFun(paramDf)
  SimData <- cbind(paramDf, SimData)
  SimData

}




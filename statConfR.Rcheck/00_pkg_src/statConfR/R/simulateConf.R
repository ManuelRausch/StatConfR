#' @title Simulate data  according to a static model of confidence

#' @description This function generates a data frame with random trials generated according to
#' the computational model of decision confidence specified in the `model` argument
#' with given parameters.
#' Simulations can be used to visualize and test qualitative model predictions
#' (e.g. using previously fitted parameters returned by \code{\link{fitConf}}).
#' See \code{\link{fitConf}} for a full mathematical description of all models
#' and their parameters.

#' @param model `character` of length 1. The generative model that should be
#'    used for simulation. Models implemented so far: 'WEV', 'SDT', 'GN', 'PDA',
#'    'IG', 'ITGc', 'ITGcm', 'logN', and 'logWEV'.
#' @param paramDf  a `data.frame`providing the number of generared trials and
#' the parameters of the chosen model. `paramDf` should contain following columns
#' (which parameters are needed depends on the specific model):
#' * \code{N} (the number of trials be simulated),
#' * \code{participant} (optional, the participant ID of each parameter set. Should be unique to each row),
#' * \code{d_1}, \code{d_2}, ... (sensitivity parameters. The number of sensitivity parameters determines the number of levels of discriminability),
#' * \code{c} (discrimination bias),
#' * \code{theta_minus.1}, \code{theta_minus.2}, ... (confidence criteria associated with the response R = -1. The function simulates one more confidence category than there are confidence criteria),
#' * \code{theta_plus.1}, \code{theta_plus.2}, ... (confidence criteria associated with the response R = 1. The function simulates one more confidence category than there are confidence criteria),
#' * \code{w} (only for models WEV and logWEV: the visibility weighting parameter, bounded between 0 and 1),
#' * \code{sigma} (only for models WEV, GN, logN, and logWEV: confidence noise, bounded between 0 and Inf),
#' * \code{m} (only for IG, ITGm, and ITGcm: metacognitive efficiency parameter, bounded between 0 and Inf),
#' * \code{b} (only for PDA: postdecisional accumulation parameter, bounded between 0 and Inf),
#' * \code{M_theta_minus.1}, \code{M_theta_minus.2}, ... (only for logN: Mean confidence criteria associated with the response R = -1),
#' * \code{M_theta_plus.1}, \code{M_theta_plus.2},... (only for logN: Mean confidence criteria associated with the response R = 1).

#' @return a dataframe with about \code{nrow(paramDf)*N} rows (see Details),
#' and the following columns:
#' - \code{participant} giving the row ID of the simulation (see Details)
#' - \code{stimulus} giving the category of the stimulus (-1 or 1)
#' - only, if more than 1 sensitivity parameter (`d1`,`d2`,...) is provided:
#' \code{diffCond} representing the difficulty condition (values correspond to
#' the levels of the sensitivity parameters, i.e. diffCond=1 represents
#' simulated trials with sensitivity `d1`)
#' - \code{response} giving the response category (-1 or 1, corresponding to the stimulus categories)
#' - \code{rating} giving the discrete confidence rating (integer, number of
#' categories depends on the number of confidence criteria provided in the parameters)
#' - \code{correct} giving the accuracy of the response (0 incorrect, 1 correct)
#' - \code{ratings} same as `rating` but as a factor

#' @details
#' The function generates about `N` trials per row with the provided parameters
#' in the data frame. The output includes a column `participant` indicating the
#' row ID of the simulated data. The values of the `participant` column may be
#' controlled by the user, by including a `participant` column in the input
#' `paramDf`. Note that the values of this column have to be unique! If no
#' `participant` column is present in the input, the row numbers will be used
#' as row IDs.
#'
#' The number of simulated trials for each row of parameters may slightly
#' deviate from the provided `N`.
#' Precisely, if there are K levels of sensitivity (i.e. there are columns
#' d1, d2, ..., dK), the function simulates `round(N/2/K)` trials per stimulus
#' identity (2 levels) and level of sensitivity (K levels).
#'
#' Simulation is performed following the generative process structure of the models.
#' See \code{fitConf} for a detailed description of the different models.

#' @author Manuel Rausch, \email{manuel.rausch@hochschule-rhein-waal.de}

#' @examples
#' # 1. define some parameters
#' paramDf <- data.frame(d_1 = 0, d_2 = 2, d_3 = 4,c = .0,
#' theta_minus.2 = -2, theta_minus.1 = -1, theta_plus.1 = 1, theta_plus.2 = 2,
#' sigma = 1/2, w = 0.5, N = 500)
#' # 2. Simulate dataset
#' SimulatedData <- simConf(model = "WEV", paramDf)

#' @importFrom plyr mdply ddply

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
  if(!"participant" %in% colnames(paramDf) ){
    paramDf$participant <- 1:nrow(paramDf)
  } else{
    if(any(duplicated(paramDf$participant))) stop("participant should be unique for each row")
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

  SimData <- plyr::ddply(paramDf, ~participant, SimFun)
  SimData
}




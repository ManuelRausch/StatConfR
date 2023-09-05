#' Fit one or more static confidence models to the data
#' @param data data.frame containing the data to fit the model to. Should contain
#' following columns: stimulus, correct, rating, and diffCond
#' @param models character vector. Models to be fit.
#' Implemented models: 'WEV', 'SDT', 'Noisy', 'PDA', and '2Chan'. Want to implement 'lognorm', 'ITGc' and 'ITGcm'
#'
#' @return data.frame with one row for each fitted model and columns for the fitted parameters and
#' as well as model fit indicates (BIC, AICc, AIC)
#' @export
#'
#' @examples
#' #To-DO

fitConfModel <- function(data ){



}

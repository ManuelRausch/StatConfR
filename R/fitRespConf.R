#' Fit a static confidence model to data
#'
#' @param data data.frame containing the data to fit the model to. Should contain
#' following columns: stimulus, correct, and rating
#' @param models character of length 1. Model to be fit.
#' Implemented models: 'WEV', 'SDT', 'Noisy', 'PDA', and '2Chan'.
#'
#' @return data.frame with one row and columns for the fitted parameters and
#' additional information (BIC, AIC,...)
#' @export
#'
#' @examples
#' #ToDo
fitRespConf <- function(data, model) {
  if (model == "WEV") {
    fitting_fct <- fitCEV
  } else if (model=="SDT") {
    fitting_fct <- fitSDT
  } else if (model=="2Chan") {
    fitting_fct <- fit2Chan
  } else if (model=="Noisy") {
    fitting_fct <- fitNoisy
  } else if (model=="PDA") {
    fitting_fct <- fitPDA
  } else stop(paste0("Model: ", model, " not implemented!\n
                                          Choose one of: 'WEV', 'SDT','2Chan','Noisy', or 'PDA'"))

  fit <- fitting_fct(data$rating, data$stimulus, data$correct, data$condition)
  return(fit)
}

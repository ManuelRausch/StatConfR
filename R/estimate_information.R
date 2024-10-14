#' Estimate Information
#'
#' Estimate the information that a classifier's responses transmit about the
#' true label. This uses a simple plug-in estimate where the conditional
#' accuracy for each response is estimated as c_i. Then, each c_i determines
#' the predictionwise (responsewise) transmitted information. The weighted
#' sum of these predictionwise transmitted information values it the overall
#' transmitted information.
#'
#' Data can be input in three ways:
#' - A data frame with variables "y" for true labels and "r" for
#'   confidence-binned responses. "y" needs to contain values -1 and +1 while
#'   r needs to be factor with levels ordered such that the first half of the
#'   factor levels are predictions for -1 and the second half for
#'   predictions +1. For example with 4 factor levels: 1: high confidence
#'   prediction for y=-1, 2: low confidence prediction for y=-1, 3: low
#'   confidence prediction for y=+1, 4: high confidence prediction for y=+1.
#' - A counts table with joint absolute frequencies. Rows correspond to true
#'   labels (stimulus categories) and columns correspond to responses.
#' - A contingency matrix with joint relative frequencies (as before but
#'   normalized to sum up to 1).
#'
#' @param x Data
#'
#'
#' @return Transmitted information

estimate_information <- function(x)
{
  UseMethod("estimate_information")
}


#'@export
estimate_information.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)

  info <- get_information(estimated_classifier)
  info
}

#'@export
estimate_information.matrix <- function(tab)
{
  estimated_classifier <- tab/sum(tab)

  info <- get_information(estimated_classifier)
  info
}

#'@export
estimate_information.table <- function(tab)
{
  estimated_classifier <- tab/sum(tab)

  info <- get_information(estimated_classifier)
  info
}

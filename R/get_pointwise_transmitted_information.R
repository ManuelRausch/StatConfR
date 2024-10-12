#' Pointwise Transmitted Information (PTI)
#' 
#' Transmitted information (often called mutual information) of a random
#' variable and response with confidence vector c.
#'
#' @param c Confidence vector, is assumed to be perfectly calibrated
#' @param p Base rate
#'
#' @return Transmitted Information between random variable and prediction
get_pointwise_transmitted_information <- pti <- function(c, prior)
{
  H(prior) - H(c)
}

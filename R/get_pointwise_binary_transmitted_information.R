#' Binary Transmitted Mutual Information (PTI)
#' 
#' Mutual information of a binary random variable (typically with equal base
#' rates, p = 0.5) and a prediction with confidence c.
#'
#' @param c Confidence in the predicted label
#' @param p Base rate
#'
#' @return Transmitted Information between binary random variable and prediction
get_pointwise_binary_transmitted_information <- pti2 <- function(c, p = 0.5) 
{ 
  H2(p) - H2(c) 
}

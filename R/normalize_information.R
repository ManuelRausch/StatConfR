#' Normalize Information
#' 
#' Given the prior and accuracy, there are limits to the transmitted
#' information. Normalize the given information by these bounds.
#' 
#' @param prior Prior of the true labels
#' @param accuracy Probability of a correct prediction (with optimal coding)
#' @param information Transmitted information
#' 
#' @return Numeric representing the normalized information (RMI)
normalize_information <- function(prior, accuracy, information)
{
  information_bounds <- get_analytic_information_bounds(prior, accuracy)
  lower_bound <- information_bounds$lowest
  upper_bound <- information_bounds$highest

  RMI <- ( information - lower_bound ) / 
         ( upper_bound - lower_bound )

  # Where bounds collapse because of accuracy edge cases, return NaN
  RMI[lower_bound == upper_bound] <- NaN
           
  RMI
}

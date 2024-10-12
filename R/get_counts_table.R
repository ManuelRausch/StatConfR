#' Get Counts Table
#' 
#' Construct a contingency table with absolute frequencies for all
#' combinations of true labels and responses.
#' 
#' @param dat Data frame with true labels representing stimulus categories
#'  (variable "y") and confidence-binned responses (variable "r")
#' 
#' @return Matrix with contingency table of absolute frequencies with true
#'  labels as rows and responses as columns
get_counts_table <- function(dat)
{
  if (!is.factor(dat$r)) stop("Variable r must be a factor with ordered levels: First half indicating y=-1 and second half indicating y=+1.")
  counts_table <- table(dat$y, dat$r)
  counts_table
}

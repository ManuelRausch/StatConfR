#' Estimate Meta-I
#'
#' Estimate Dayan (2023)'s information-theoretic measure of metacognitive accuracy meta-I.
#' Meta-I is the transmitted information minus the minimal information given the accuracy of the classifier,
#' $$meta-I = I(Y; \hat{Y}, C) - I(Y; \hat{Y})$$
#'
#' This is equivalent to Dayan's formulation where meta-I is the information
#' that confidences transmit about the correctness of a response,
#' $$meta-I = I(Y = \hat{Y}; C).$$
#'
#' Data can be input in three ways:
#' - A data frame with variables "y" for true labels and "r" for
#'   confidence-binned responses. "y" needs to contain values -1 and +1 while
#'   r needs to be a factor with ordered levels such that the first half of
#'   the levels correspond to predictions for y=-1 and the second half to
#'   predictions for y=+1.
#' - A counts table with joint absolute frequencies. Rows correspond to true
#'   labels (stimulus categories) and columns correspond to responses.
#' - A contingency matrix with joint relative frequencies (as before but
#'   normalized to sum up to 1).
#'
#' @return Meta-I value
#' @author Sascha Meyen
#' @name estimate_meta_I (expressed in bits, i.e. log base is 2)
#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
#' @export
estimate_meta_I <- function(x)
{
  UseMethod("estimate_meta_I")
}

#' @exportS3method
estimate_meta_I.matrix <- function(estimated_classifier)
{
  estimated_classifier <- estimated_classifier/sum(estimated_classifier)
  p                    <- rowSums(estimated_classifier) # Prior

  info       <- get_information(estimated_classifier)
  a          <- get_accuracy(estimated_classifier)
  info_lower <- get_lower_info_for_one(prior    = p,
                                       accuracy = a)

  meta_I <- info - info_lower
  meta_I
}


#' @exportS3method
estimate_meta_I.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)

  meta_I <- estimate_meta_I.matrix(estimated_classifier)
  meta_I
}


#' @exportS3method
estimate_meta_I.table <- function(counts_table)
{
  estimated_classifier <- counts_table/sum(counts_table)

  meta_I <- estimate_meta_I.matrix(estimated_classifier)
  meta_I
}

#' Estimate Meta-I
#'
#' @description Estimate meta-I, an information-theoretic measure of metacognitive sensitivity
#' proposed by Dayan (2023).
#'
#' @params x
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
#' @details
#' Meta-I is the mutual information between the confidence and accuracy and is calculated as
#' the transmitted information minus the minimal information given the accuracy,
#' /deqn{meta-I = I(Y; \hat{Y}, C) - I(Y; \hat{Y})}. This is equivalent to Dayan's formulation where meta-I is the information
#' that confidences transmit about the correctness of a response, /deqn{meta-I = I(Y = \hat{Y}; C).}

#' @return Meta-I value (expressed in bits, i.e. log base is 2)

#' @author Sascha Meyen

#' @name estimate_meta_I
#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
#' @export
estimate_meta_I <- function(x)
{
  UseMethod("estimate_meta_I")
}

#' @keywords internal
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

#' @keywords internal
estimate_meta_I.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)

  meta_I <- estimate_meta_I.matrix(estimated_classifier)
  meta_I
}

#' @keywords internal
estimate_meta_I.table <- function(counts_table)
{
  estimated_classifier <- counts_table/sum(counts_table)

  meta_I <- estimate_meta_I.matrix(estimated_classifier)
  meta_I
}

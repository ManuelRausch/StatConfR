#' Estimate Meta Ir2
#'
#' Estimate Dayan (2023)'s meta-I^r_1. This is meta-I (see
#' estimate_meta_I()) normalized by the remaining uncertainty
#' about whether the predictions are correct or not, H(Y = \hat{Y}),
#'
#'   meta-Ir2 = meta-I / H(Y = \hat{Y}).
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
#' @return Meta-I^r_2 value
#' @author Sascha Meyen
#' @name estimate_meta_Ir2
#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
#' @export
estimate_meta_Ir2 <- function(x)
{
  UseMethod("estimate_meta_Ir2")
}

estimate_meta_Ir2.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)
  meta_I <- estimate_meta_Ir2.matrix(estimated_classifier)
  meta_I
}

estimate_meta_Ir2.matrix <- function(estimated_classifier)
{
  estimated_classifier <- estimated_classifier/sum(estimated_classifier)

  meta_I <- estimate_meta_I(estimated_classifier)

  a          <- get_accuracy(estimated_classifier)
  H_accuracy <- H2(a)

  meta_I_r2 <- meta_I / H_accuracy

  meta_I_r2
}

estimate_meta_Ir2.table <- function(counts_table)
{
  estimated_classifier <- counts_table/sum(counts_table)
  meta_Ir2 <- estimate_meta_Ir2.matrix(estimated_classifier)
  meta_Ir2
}

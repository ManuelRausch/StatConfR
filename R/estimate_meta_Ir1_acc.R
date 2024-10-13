#' Estimate Meta Ir1 via Accuracy
#'
#' Estimate Dayan (2023)'s information-theoretic measure of metacognitive sensitivity, meta-I^r_1.
#' This is meta-I (see estimate_meta_I()) normalized by the meta-I that would be
#' produced by normal noise with the estimated sensitivity,
#'
#'   $$meta-Ir1 = meta-I / meta-I(d')$$.
#'
#' In contrast to estimate_meta_Ir1, the sensitivity is here computed by
#' determining the accuracy and converting it to a sensitivity assuming an
#' unbiased observer, d' = 2*qnorm(accuracy).
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
#' @return Meta-I^r_1 value
#'
#' @export
estimate_meta_Ir1_acc <- function(x, ...)
{
  UseMethod("estimate_meta_Ir1_acc")
}

#' @keywords internal
estimate_meta_Ir1_acc.matrix <- function(counts_table)
{
  estimated_classifier <- counts_table/sum(counts_table)

  # (Unnormalized) meta-I
  meta_I        <- estimate_meta_I(estimated_classifier)

  # meta-I of a normal noise classifier
  a             <- estimate_accuracy(estimated_classifier)
  d             <- 2*qnorm(a)
  info_normal   <- get_normal_noise_information(sensitivities = d)$info
  p             <- rowSums(estimated_classifier)/sum(estimated_classifier)
  info_lower    <- get_lower_info_for_one(p, a)
  meta_I_normal <- info_normal - info_lower

  # Normalization
  meta_I_r1 <- meta_I / meta_I_normal

  meta_I_r1
}

#' @keywords internal
estimate_meta_Ir1_acc.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)
  meta_I <- estimate_meta_Ir1_acc.matrix(estimated_classifier)
  meta_I
}
#' @keywords internal
estimate_meta_Ir1_acc.table <- function(counts_table)
{
  meta_Ir1 <- estimate_meta_Ir1_acc.matrix(counts_table)
  meta_Ir1
}

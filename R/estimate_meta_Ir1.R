#' Estimate Meta Ir1
#'
#' Estimate Dayan (2023)'s information-theoretic measure of metacognitive accuracy meta-I^r_1.
#' This is meta-I (see \code{\link{ estimate_meta_I}})
#' normalized by meta-I that would be expected from a signal detection rating model
#'
#'   meta-Ir1 = meta-I / meta-I(d').
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
#' @author Sascha Meyen
#' @name estimate_meta_Ir1
#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
#' @export
estimate_meta_Ir1 <- function(x, ...)
{
  UseMethod("estimate_meta_Ir1")
}

estimate_meta_Ir1.matrix <- function(counts_table)
{
  estimated_classifier <- counts_table/sum(counts_table)

  # (Unnormalized) meta-I
  meta_I        <- estimate_meta_I(estimated_classifier)

  # meta-I of a normal noise classifier
  p             <- rowSums(estimated_classifier)
  d             <- estimate_sensitivity(counts_table)
  a             <- pnorm(abs(d)/2)
  info_lower    <- get_lower_info_for_one(p, a)
  info_normal   <- get_normal_noise_information(sensitivities = d)$info
  meta_I_normal <- info_normal - info_lower

  # Normalization
  meta_I_r1 <- meta_I / meta_I_normal

  meta_I_r1
}

estimate_meta_Ir1.data.frame <- function(msd)
{
  counts_table <- get_counts_table(msd)
  meta_I <- estimate_meta_Ir1.matrix(counts_table)
  meta_I
}

estimate_meta_Ir1.table <- function(counts_table)
{
  meta_Ir1 <- estimate_meta_Ir1.matrix(counts_table)
  meta_Ir1
}

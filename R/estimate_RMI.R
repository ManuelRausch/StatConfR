#' Estimate Relative Meta-Information
#'
#' Estimate relative meta-information, RMI. This is meta-I (see
#' estimate_meta_I()) normalized by the possible range of meta-I values
#' for the estimated accuracy,
#'
#'   RMI = meta-I / max meta-I.
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
#' @param x Data
#'
#' @return Relative meta-information value
#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
#' @author Sascha Meyen
#' @name estimate_RMI
#' @export
estimate_RMI <- function(x)
{
  UseMethod("estimate_RMI")
}

estimate_RMI.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)
  estimated_RMI <- estimate_RMI.matrix(estimated_classifier)
  estimated_RMI
}

estimate_RMI.matrix <- function(estimated_classifier)
{
  estimated_classifier <- estimated_classifier/sum(estimated_classifier)
  estimated_RMI <- get_RMI(estimated_classifier)
  estimated_RMI
}

estimate_RMI.table <- function(tab)
{
  estimated_classifier <- tab/sum(tab)
  estimated_RMI <- get_RMI(estimated_classifier)
  estimated_RMI
}

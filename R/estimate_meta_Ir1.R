#' @title Estimate meta-\eqn{I_{1}^{r}}

#' @description
#' `estimate_meta_Ir1` estimates meta-\eqn{I_{1}^{r}}, an information-theoretic measure
#'  of metacognitive efficiency proposed by Dayan (2023), using the sensitivity index d'.

#' @details Meta-\eqn{I_{1}^{r}} is meta-I (see \code{\link{estimate_meta_I}})
#' normalized by the value of meta-I expected assuming a signal detection model (Green & Swets, 1966)
#' with Gaussian noise, based on calculating the sensitivity index d':
#' \deqn{meta-I_{1}^{r} = meta-I / meta-I(d')}
#' It should be noted that Dayan (2023) pointed out that a liberal or conservative use of the confidence levels
#' will affected the mutual information and thus influence meta-\eqn{I_{1}^{r}}.
#'
#' @param x Three different types of inputs are accepted:
#' - A `data.frame` with variables "y" for true labels and "r" for
#'   confidence-binned responses. "y" needs to contain values -1 and +1 while
#'   r needs to be a factor with ordered levels such that the first half of
#'   the levels correspond to predictions for y = -1 and the second half to
#'   predictions for y = +1.
#' - A counts `table` with joint absolute frequencies. Rows correspond to true
#'   labels (stimulus categories) and columns correspond to responses.
#' - A contingency `matrix` with joint relative frequencies (as before but
#'   normalized to sum up to 1).

#' @return meta-\eqn{I_{1}^{r}} value (expressed in bits, i.e. log base is 2)

#' @examples
#' # 1. prepare counts table for one subject
#' OneSbj <- subset(MaskOri, participant == 1)
#' y <- ifelse(OneSbj$stimulus == 0, -1, 1)
#' r <- factor(ifelse(OneSbj$response == 0, -1, 1) * as.numeric(OneSbj$rating))
#' counts <- table(y, r)
#'
#' # 2. calculate meta-I
#' meta_Ir1 <- estimate_meta_Ir1(counts)

#' @author Sascha Meyen, \email{saschameyen@gmail.com}

#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
#' @references Green, D. M., & Swets, J. A. (1966). Signal detection theory and psychophysics. Wiley.\cr
#'
#' @export
estimate_meta_Ir1 <- function(x)
{
  UseMethod("estimate_meta_Ir1")
}
#' @export
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

#' @export
estimate_meta_Ir1.data.frame <- function(msd)
{
  counts_table <- get_counts_table(msd)
  meta_I <- estimate_meta_Ir1.matrix(counts_table)
  meta_I
}

#' @export
estimate_meta_Ir1.table <- function(counts_table)
{
  meta_Ir1 <- estimate_meta_Ir1.matrix(counts_table)
  meta_Ir1
}

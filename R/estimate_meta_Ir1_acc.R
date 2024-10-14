#' @title Estimate Estimate meta-\eqn{I_{1}^{r}} via Accuracy

#' @description `estimate_meta_Ir1_acc` estimates meta-\eqn{I_{1}^{r}},
#' an information-theoretic measure of metacognitive efficiency proposed by Dayan (2023),
#' using the proportion of correct responses.

#' @details Meta-\eqn{I_{1}^{r}} is meta-I (see \code{\link{estimate_meta_I}})
#' normalized by the value of meta-I expected assuming a signal detection model (Green & Swets, 1966)
#' with Gaussian noise as well as an unbiased decision criterion:
#' \deqn{meta-I_{1}^{r} = meta-I / meta-I(d')}
#' In contrast to \code{\link{estimate_meta_Ir1}}, `estimate_meta_Ir1_acc` computes the  sensitivity by
#' determining the proportion of correct responses and then converting it to a sensitivity assuming an
#' unbiased observer, \eqn{d' = 2 \times \Phi^{-1}(P_{correct})}..
#' It should also be noted that Dayan (2023) pointed out that a liberal or conservative use of the confidence levels
#' will affected the mutual information and thus influence meta-\eqn{I_{1}^{r}}.

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
#' meta_Ir1 <- estimate_meta_Ir1_acc(counts)

#' @author Sascha Meyen, \email{saschameyen@gmail.com}

#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
#' @references Green, D. M., & Swets, J. A. (1966). Signal detection theory and psychophysics. Wiley.\cr

#' @export
estimate_meta_Ir1_acc <- function(x, ...)
{
  UseMethod("estimate_meta_Ir1_acc")
}

#'@export
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

#'@export
estimate_meta_Ir1_acc.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)
  meta_I <- estimate_meta_Ir1_acc.matrix(estimated_classifier)
  meta_I
}
#'@export
estimate_meta_Ir1_acc.table <- function(counts_table)
{
  meta_Ir1 <- estimate_meta_Ir1_acc.matrix(counts_table)
  meta_Ir1
}

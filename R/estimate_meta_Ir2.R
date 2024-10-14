#' @title Estimate meta-\eqn{I_{2}^{r}}

#' @description
#' `estimate_meta_Ir2` estimates meta-\eqn{I_{2}^{r}}, an information-theoretic measure
#'  of metacognitive efficiency proposed by Dayan (2023).

#' @details
#' Meta-\eqn{I_{2}^{r}} is meta-I (see \code{\link{estimate_meta_I}}) normalized
#' by its theoretical upper bound, which is the information entropy of accuracy,  \eqn{H(Y = \hat{Y})}:
#' \deqn{meta-I_{2}^{r} = meta-I / H(Y = \hat{Y})}
#' It should be noted that Dayan (2023) pointed out that a liberal or conservative use of the confidence levels
#' will affected the mutual information and thus influence meta-\eqn{I_{2}^{r}}.

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

#' @return meta-\eqn{I_{2}^{r}} value (expressed in bits, i.e. log base is 2)

#' @examples
#' # 1. prepare counts table for one subject
#' OneSbj <- subset(MaskOri, participant == 1)
#' y <- ifelse(OneSbj$stimulus == 0, -1, 1)
#' r <- factor(ifelse(OneSbj$response == 0, -1, 1) * as.numeric(OneSbj$rating))
#' counts <- table(y, r)
#'
#' # 2. calculate meta-Ir2
#' meta_Ir2 <- estimate_meta_Ir2(counts)

#' @author Sascha Meyen, \email{saschameyen@gmail.com}
#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>

#' @export
estimate_meta_Ir2 <- function(x){
  UseMethod("estimate_meta_Ir2")
}

#' @export
estimate_meta_Ir2.data.frame <- function(x){
  estimated_classifier <- estimate_classifier(x)
  meta_I <- estimate_meta_Ir2.matrix(estimated_classifier)
  meta_I
}

#' @export
estimate_meta_Ir2.matrix <- function(x){
  estimated_classifier <- x/sum(x)

  meta_I <- estimate_meta_I(estimated_classifier)

  a          <- get_accuracy(estimated_classifier)
  H_accuracy <- H2(a)

  meta_I_r2 <- meta_I / H_accuracy

  meta_I_r2
}

#' @export
estimate_meta_Ir2.table <- function(x){
  estimated_classifier <- x/sum(x)
  meta_Ir2 <- estimate_meta_Ir2.matrix(estimated_classifier)
  meta_Ir2
}

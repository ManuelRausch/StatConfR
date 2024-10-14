#' @title Estimate relative Meta-Information

#' @description `estimate_RMI` estimates the RMI index of relative meta-information.

#' @details
#' RMI is meta-I (see \code{\link{estimate_meta_I}}) normalized by the maximum
#' possible range of meta-I possible for the estmiated level of accuracy:
#' \deqn{RMI = meta-I / max(meta-I)}
#' It should be noted that Dayan (2023) pointed out that a liberal or conservative use of the confidence levels
#' will affected the mutual information and thus influence RMI.

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

#' @return Relative meta-information value, expressed in bits (i.e. log base is 2).
#'
#' @examples
#' # 1. prepare counts table for one subject
#' OneSbj <- subset(MaskOri, participant == 1)
#' y <- ifelse(OneSbj$stimulus == 0, -1, 1)
#' r <- factor(ifelse(OneSbj$response == 0, -1, 1) * as.numeric(OneSbj$rating))
#' counts <- table(y, r)
#'
#' # 2. calculate meta-Ir2
#' RMI <- estimate_RMI(counts)

#' @author Sascha Meyen, \email{saschameyen@gmail.com}
#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>

#' @export
estimate_RMI <- function(x){
  UseMethod("estimate_RMI")
}

#' @export
estimate_RMI.data.frame <- function(x){
  estimated_classifier <- estimate_classifier(x)
  estimated_RMI <- estimate_RMI.matrix(estimated_classifier)
  estimated_RMI
}

#' @export
estimate_RMI.matrix <- function(x){
  estimated_classifier <- x/sum(x)
  estimated_RMI <- get_RMI(estimated_classifier)
  estimated_RMI
}

#' @export
estimate_RMI.table <- function(x){
  estimated_classifier <- x/sum(x)
  estimated_RMI <- get_RMI(estimated_classifier)
  estimated_RMI
}

#' @title Estimate Meta-I

#' @description `estimate_meta_I` estimates meta-I, an information-theoretic measure
#'  of metacognitive sensitivity proposed by Dayan (2023).

#' @details
#' Meta-I is defined as the mutual information between the confidence and accuracy and is calculated as
#' the transmitted information minus the minimal information given the accuracy,
#' \deqn{meta-I = I(Y; \hat{Y}, C) - I(Y; \hat{Y})}.
#' This is equivalent to Dayan's formulation where meta-I is the information
#' that confidences transmit about the correctness of a response,
#' \deqn{meta-I = I(Y = \hat{Y}; C).}
#' It should be noted that Dayan (2023) pointed out that a liberal or conservative use of the confidence levels
#' will affected th mutual information and thus influence meta-I.
#'
#' @param x Three different types of inputs are accepted:
#' - A `data.frame` with variables "y" for true labels and "r" for
#'   confidence-binned responses. "y" needs to contain values -1 and +1 while
#'   r needs to be a factor with ordered levels such that the first half of
#'   the levels correspond to predictions for y=-1 and the second half to
#'   predictions for y=+1.
#' - A counts `table` with joint absolute frequencies. Rows correspond to true
#'   labels (stimulus categories) and columns correspond to responses.
#' - A contingency `matrix` with joint relative frequencies (as before but
#'   normalized to sum up to 1).

#' @return Meta-I value (expressed in bits, i.e. log base is 2)

#' @examples
#' # 1. prepare counts table for one subject
#' OneSbj <- subset(MaskOri, participant == 1)
#' y <- ifelse(OneSbj$stimulus == 0, -1, 1)
#' r <- factor(ifelse(OneSbj$response == 0, -1, 1) * as.numeric(OneSbj$rating))
#' counts <- table(y, r)
#'
#' # 2. calculate meta-I
#' metaI <- estimate_meta_I(counts)
#'
#' @author Sascha Meyen, \email{saschameyen@gmail.com}

#' @references Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
# to do: add Saschas Paper once it is available ;-)

#' @export
estimate_meta_I <- function(x){
  UseMethod("estimate_meta_I")
}

#'@export
estimate_meta_I.matrix <- function(x){
  estimated_classifier <- x/sum(x)
  p                    <- rowSums(x) # Prior

  info       <- get_information(x)
  a          <- get_accuracy(x)
  info_lower <- get_lower_info_for_one(prior    = p,
                                       accuracy = a)
  meta_I <- info - info_lower
  meta_I
}

#' @export
estimate_meta_I.data.frame <- function(x){
  estimated_classifier <- estimate_classifier(x)
  meta_I <- estimate_meta_I.matrix(estimated_classifier)
  meta_I
}

#' @export
estimate_meta_I.table <- function(x){
  estimated_classifier <- x/sum(x)
  meta_I <- estimate_meta_I.matrix(estimated_classifier)
  meta_I
}

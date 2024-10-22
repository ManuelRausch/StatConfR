#' @title Estimate Metacognitive Information Theory Measures

#' @description `estimateMetaI` estimates meta-\eqn{I}, an information-theoretic
#'  measure of metacognitive sensitivity proposed by Dayan (2023), as well as
#'  derived measures like it including meta-\eqn{I_{1}^{r}} and Meta-\eqn{I_{2}^{r}}.
#'  These are different normalizations of meta-\eqn{I}:
#' - meta-\eqn{I_{1}^{r}} normalizes by the meta-\eqn{I} that would be
#'   expected from an underlying normal distribution with the same
#'   sensitivity.
#' - meta-\eqn{I_{1}^{r\prime}} is a variant not discussed by Dayan
#'   (2023) which normalizes by the meta-\eqn{I} that would be expected from
#'   an underlying normal distribution with the same accuracy.(This is
#'   similar to the sensitivity approach but without considering variable
#'   thresholds.)
#' - meta-\eqn{I_{2}^{r}} normalizes by the maximum amount of meta-\eqn{I}
#'   which would be reached if all uncertainty about the stimulus was removed.
#'  Additionally, one more measure that was not discussed by Dayan (2023) is
#'  provided.
#' - \eqn{RMI} normalizes meta-\eqn{I} by the range of its possible
#'    values and therefore scales between 0 and 1.
#'  All measures are also calculated with a bias-reduced variant for which the
#'  observed frequencies are taken as underlying probability distribution to
#'  estimate the sampling bias. The estimated bias is then subtraced from the
#'  initial measures. This approach uses Monte-Carlo simulations and is
#'  therefore not deterministic (values can vary from one evaluation of the
#'  function to the next). However, it is a simple way to reduce the bias
#'  inherent in these measures.

#' @details
#' Meta-\deqn{I} is defined as the mutual information between the confidence and
#' accuracy and is calculated as the transmitted information minus the
#' minimal information given the accuracy,
#' meta-\deqn{I = I(Y; \hat{Y}, C) - I(Y; \hat{Y})}.
#' This is equivalent to Dayan's formulation where meta-I is the information
#' that confidences transmit about the correctness of a response,
#' meta-\deqn{I = I(Y = \hat{Y}; C).}
#' The other measures here are different normalizations of this Meta-\deqn
#' {I}. It should be noted that Dayan (2023) pointed out that a liberal or
#' conservative use of the confidence levels will affected th mutual
#' information and thus influence meta-I.
#'
#' @param data  a `data.frame` where each row is one trial, containing following
#' variables:
#' * \code{participant} (some group ID, most often a participant identifier;
#'    the meta-I measures are estimated for each subset of `data`
#'    determined by the different values of this column),
#' * \code{stimulus} (stimulus category in a binary choice task,
#'    should be a factor with two levels, otherwise it will be transformed to
#'    a factor with a warning),
#' * \code{rating} (discrete confidence judgments, should be a factor with levels
#'    ordered from lowest confidence to highest confidence;
#'    otherwise will be transformed to factor with a warning),
#' * \code{correct} (encoding whether the response was correct; should  be 0 for
#'    incorrect responses and 1 for correct responses)
#' @param bias_reduction `logical`. Whether to apply the bias reduction or
#'  not. If runtime is to long, consider setting this to FALSE
#'  (default: TRUE).

#' @return Estimated Meta-\eqn{I} value (expressed in bits, i.e. log base is
#'  2, or unitless for the normalized measures)

#' @examples
#' # 1. Select two subjects from the masked orientation discrimination experiment
#' data <- subset(MaskOri, participant %in% c(1:2))
#' head(data)
#'
#' # 2. Calculate meta-I measures with bias reduction (this may take 10 s per subject)
#' \donttest{
#' metaIMeasures <- estimateMetaI(data)
#' metaIMeasures
#' }
#'
#' # 3. Calculate meta-I measures for all participants without bias reduction (much faster)
#' metaIMeasures <- estimateMetaI(MaskOri, bias_reduction = FALSE)
#' metaIMeasures

#' @author Sascha Meyen, \email{saschameyen@gmail.com}

#' @references Dayan, P. (2023). Metacognitive Information Theory.
#'  Open Mind, 7, 392–411. <https://doi.org/10.1162/opmi_a_00091>
# to do: add Saschas Paper once it is available ;-)

#' @export
estimateMetaI <- function(data, bias_reduction = TRUE) {

  # Enforce relevant data variables to be factors
  no_participants_variable <- FALSE
  if (is.null(data$participant)) {
    data$participant <- 1
    no_participants_variable <- TRUE
  }

  if (!all(c("stimulus", "rating", "correct") %in% names(data) )) {
    data$stimulus <- factor(data$stimulus)
    warning("stimulus is transformed to a factor!")
  }

  if (!is.factor(data$stimulus)) {
    data$stimulus <- factor(data$stimulus)
    warning("stimulus is transformed to a factor!")
  }

  if(length(unique(data$stimulus)) != 2) {
    stop("There must be exacltly two different possible values of stimulus")
  }

  if (!is.factor(data$rating)) {
    data$rating <- factor(data$rating)
    warning("rating is transformed to a factor!")
  }

  if(!all(data$correct %in% c(0,1))) {
    stop("correct should be 1 or 0")
  }

  # Estimate information-theoretic measures of metacognition
  res <- data.frame()
  for (participant in unique(data$participant))
  {
    s <- participant == data$participant
    part_data <- data[s, ]

    estimated_classifier <- estimate_classifier(part_data)
    number_of_stimuli    <- table(part_data$stimulus)

    re <- data.frame(participant           = participant                                ,
                     meta_I                = estimate_meta_I(estimated_classifier)      ,
                     meta_Ir1              = estimate_meta_Ir1(estimated_classifier)    ,
                     meta_Ir1_acc          = estimate_meta_Ir1_acc(estimated_classifier),
                     meta_Ir2              = estimate_meta_Ir2(estimated_classifier)    ,
                     RMI                   = estimate_RMI(estimated_classifier)         )
    if (bias_reduction)
    {
      rb <- data.frame(meta_I_debiased       = get_bias_reduced_meta_measure(estimated_classifier,
                                                                             number_of_stimuli   ,
                                                                             estimate_meta_I     ) ,
                       meta_Ir1_debiased     = get_bias_reduced_meta_measure(estimated_classifier,
                                                                             number_of_stimuli   ,
                                                                             estimate_meta_Ir1   ) ,
                       meta_Ir1_acc_debiased = get_bias_reduced_meta_measure(estimated_classifier ,
                                                                             number_of_stimuli    ,
                                                                             estimate_meta_Ir1_acc),
                       meta_Ir2_debiased     = get_bias_reduced_meta_measure(estimated_classifier,
                                                                             number_of_stimuli   ,
                                                                             estimate_meta_Ir2   ) ,
                       RMI_debiased          = get_bias_reduced_meta_measure(estimated_classifier,
                                                                             number_of_stimuli   ,
                                                                             estimate_RMI        ) )
      re <- cbind(re, rb)
    }
    res <- rbind(res, re)
  }

  # Drop unnecessary variables
  if (no_participants_variable) res$participant <- NULL

  res
}

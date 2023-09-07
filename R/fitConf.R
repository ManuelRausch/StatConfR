#' Fit a static confidence model to data
#'
#' This function fits one model of decision confidence to a data set.
#' It calls a corresponding fitting function for the selected model. For full details,
#' see Rausch et al. (2018), Rausch et al. (2020) or Rausch et al. (2023).
#'
#' @param data  a `data.frame` where each row is one trial, containing following
#' variables:
#' * \code{condition} (optional; for different levels of stimulus quality,
#'    should be a factor with levels ordered from hardest to easiest),
#' * \code{rating} (discrete confidence judgments, should be given as factor;
#'    otherwise will be transformed to factor with a warning),
#' * \code{stimulus} (encoding the stimulus category in a binary choice task,
#'    should be a factor with two levels, otherwise it will be transformed to
#'    a factor with a warning),
#' * \code{correct} (encoding whether the response was correct; values in 0, 1)
#' @param model `character` of length 1.
#' Implemented models: 'WEV', 'SDT', 'Noisy', 'PDA', 'IG', 'ITGc' and *ITGcm'
#' Alternatively, if `model="all"` (default), all implemented models will be fit.
#' @param var `character`. One of "constant" (default), "increasing", or "free" (will be
#' implemented), indicating how noise variances should be treated across conditions.
#' See Details for more information. Applies only to the models "SDT" and "WEV"!
#'
#' @return Gives data frame with one row and columns for the fitted parameters of the
#' selected model as well as additional information about the fit (`negLogLik` (for final parameters),
#' `k` (number of parameters), `N` (number of data rows), `BIC`, `AICc` and `AIC`)
#'
#'#' @details The fitting involves a first grid search through an initial grid. Then the best 25 (\code{nAttempts})
#' parameter sets are chosen for an optimization, which is done with the Nelder-Mead algorithm implemented in \code{\link[stats]{optim}}.
#'
#' ## Mathematical description of models
#' This section contains a detailed mathematical description of all models
#' implemented in the package.
#'
#' The computational models are all based on signal detection theory. Assume that
#' there are \eqn{k} different levels of difficulty manipulated in the
#' (and the levels are given by the `condition` column in the `data`) and that
#' the `stimulus` column indicated the identity of the true stimulus \eqn{S}
#' being either -1 or 1. Then, for each level of difficulty , a value for the
#' sensitivity \eqn{d_k} is fit. The models assume that the stimulus
#' generates normally distributed sensory evidence \eqn{x} with mean \eqn{Sd_k/2}
#' and variance \eqn{\sigma_k} (see below). The sensory evidence \eqn{x}
#' is compared to a decision threshold \eqn{\theta} to generate a choice response
#' \eqn{R}, which is 1, if \eqn{x} exceeds \eqn{\theta} and -1 else. (In the
#' output of the functions this will be A and B respectively.)
#' To generate confidence, the confidence variable \eqn{y} is compared to another
#' set of thresholds \eqn{c_{D,i}, D=A, B,  i=1,...,L-1}, depending on the
#' initial choice \eqn{D} to produce a \eqn{L}-step discrete confidence response.
#' The number of thresholds will be inferred by the number of steps in the
#' `rating` column of `data`.
#' The parameters common to all models are thus:
#' - sensitivity parameters \eqn{d_1},...,\eqn{d_k} (\eqn{k}: number of difficulty levels)
#' - choice threshold \eqn{\theta}
#' - confidence threshold \eqn{c_{A,1}},...\eqn{c_{A,L-1}},\eqn{c_{B,1}},...
#' \eqn{c_{B,L-1}} (\eqn{L}: number of steps for confidence ratings)
#'
#' How the confidence variable \eqn{y} is computed
#' varies from model to model. Following models are implemented:
#'
#' ### \strong{signal-detection theory (SDT)}
#' According to the signal-detection theory for confidence, the same sensory
#' evidence used to generate the response is used to generate confidence, i.e.
#' \eqn{y=x} and the confidence thresholds span from the left and
#' right side of the decision threshold \eqn{\theta}.
#'
#' ### \strong{Noisy signal-detection theory (Noisy)}
#' According to the noisy signal-detection theory, \eqn{y} is subject to
#' additional noise and assumed to be normally distributed around the initial
#' evidence value \eqn{x} with some standard deviation \eqn{\sigma}.
#' \eqn{\sigma} is an additional parameter fit to the model.
#'
#' ### \strong{weighted evidence and visibility model (WEV)}
#' WEV assumes that the observer combines evidence about choice-relevant features
#' of the stimulus with the strength of evidence about choice-irrelevant features
#' to generate confidence. Thus, the WEV model assumes that \eqn{y} is normally
#' distributed with a mean of \eqn{(1-w)x+wdR} and standard deviation \eqn{\sigma}.
#' The standard deviation quantifies the amount of unsystematic variability
#' contributing to confidence judgments but not to the identification judgments.
#' The parameter \eqn{w} represents the weight that is put on the choice-irrelevant
#' features in the confidence judgment. \eqn{w} and \eqn{\sigma} are fitted in
#' addition to the common parameters.
#'
#' ### \strong{Post-decisional accumulation model (PDA)}
#' PDA incorporates the idea of ongoing information accumulation after the
#' initial choice. The parameter \eqn{a} indicates the time of additional
#' accumulation. The confidence variable is normally distributed with mean
#' \eqn{x+Sda} and variance \eqn{a}.
#' For this model the parameter \eqn{a} is fitted in addition to the common
#' parameters.
#'
#' ### \strong{Two-channel model (2Chan)}
#' According to the two-channel model, \eqn{y} is sampled completely independent
#' from \eqn{x}. It is normally distributed with a mean of \eqn{ad} and variance
#' of 1 (again as it would scale with \eqn{a}). The additional parameter \eqn{a}
#' represents the amount of information available for the metacognitive judgment
#' relative to the type-I choice and can be smaller as well as greater than 1.
#'
#' ### `var` argument
#' The `var`argument will only apply to SDT and WEV models.
#' Depending on the `var` argument, the variance of the sensory evidence
#' \eqn{\sigma_k^2} will be treated as:
#' - `constant` across conditions. In this case \eqn{\sigma_k^2 = 1} for all \eqn{k} since
#' it would scale with the sensitivity parameters
#' - `increasing` across conditions. In this case
#' \eqn{\sigma_k^2 = 1+a(\frac{d_k}{2})^2} would increase with
#' sensitivity across conditions. \eqn{a} is a parameter that will be fit additionally.
#'
#' @md
#'
#' @author Sebastian Hellmann, \email{sebastian.hellmann@@ku.de}
#'
#' @name fitConfModel
#' @importFrom stats dnorm pnorm optim integrate
#'
#' @references Rausch, M., Hellmann, S. & Zehetleitner, M. (2018). Confidence in
#'  masked orientation judgments is informed by both evidence and visibility.
#'  \emph{Atten Percept Psychophys} 80, 134–154. doi: 10.3758/s13414-017-1431-5
#'
#' @examples
#' # 1. Generate data from an artificial participant
#' ###      d1,  d2, d3, cA3, cA2, cA1,theta, cB1, cB2, cB3)
#' p1 <- c(0.2, 0.5, 1.5,-1.0,-0.5,-0.2,  0.1, 0.3, 0.6, 0.9)
#' D <- p1[1:3]
#' thresholds <- c(-Inf, p1[4:10], Inf)
#' data <- expand.grid(n = 1:30, stimulus = c(-1, 1), condition=c(1,2, 3))
#' data$x <- rnorm(nrow(data), mean=data$stimulus*D[data$condition])
#' data$response1 <- as.numeric(cut(data$x,breaks=thresholds, include.lowest = TRUE))
#' data$response <- ifelse(data$response1<=4, -1, 1)
#' data$correct <- as.numeric(data$response==data$stimulus)
#' data$rating <- ifelse(data$response==-1, 5-data$response1, data$response1-4)
#' data$stimulus <- as.factor(data$stimulus)
#' data$condition <- as.factor(data$condition)
#' table(data[data$correct==1,c("condition")])/60
#' head(data)
#'
#'
#' # 2. Use fitting function
#' \dontrun{
#'   # Fitting takes some time to run:
#'   fitConf(data, model=c("WEV"))
#' }
#'
#'

#' @export
fitConf <- function(data, model#, var="constant"
) {
  if (is.null(data$condition)) data$condition <- 1
  if (!is.factor(data$condition)) {
    data$condition <- factor(data$condition)
    warning("condition transformed to a factor!")
  }
  if (!is.factor(data$stimulus)) {
    data$stimulus <- factor(data$stimulus)
    warning("stimulus transformed to a factor!")
  }
  if (!is.factor(data$rating)) {
    data$rating <- factor(data$rating)
    warning("rating  transformed to a factor!")
  }

  if (model == "WEV") {
    fitting_fct <- fitCEV
    #if (var=="increasing") fitting_fct <- fitCEVvarS
  } else if (model=="SDT") {
    fitting_fct <- fitSDT
    #if (var=="increasing") fitting_fct <- fitSDTvarS
  } else if (model=="IG") {
    fitting_fct <- fit2Chan
  } else if (model=="ITGc") {
    fitting_fct <- fitITGc
  } else if (model=="ITGcm") {
    fitting_fct <- fitITGcm
  } else if (model=="Noisy") {
    fitting_fct <- fitNoisy
  } else if (model=="PDA") {
    fitting_fct <- fitPDA
  } else stop(paste0("Model: ", model, " not implemented!\nChoose one of: 'WEV', 'SDT','IG', 'ITGc', 'ITGcm,'Noisy', or 'PDA'"))

  fit <- fitting_fct(data$rating, data$stimulus, data$correct, data$condition)
  return(fit)
}

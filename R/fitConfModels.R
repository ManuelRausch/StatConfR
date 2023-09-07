#' Fit several static confidence models to multiple participants
#'
#' This function is a wrapper of the function \code{\link{fitConf}} (see
#' there for more information). It calls the function for every possible combination
#' of model in the `model` argument and participant in the \code{data}, respectively.
#' See the Details for more information about the parameters. All models are
#' described in full detail in Rausch et al. (2018).
#'
#' @param data  a `data.frame` where each row is one trial, containing following
#' variables:
#' * \code{condition} (not necessary; for different levels of stimulus quality,
#'    should be a factor, otherwise it will be transformed to a factor with a
#'    waring),
#' * \code{rating} (discrete confidence judgments, should be given as factor;
#'    otherwise will be transformed to factor with a warning),
#' * \code{stimulus} (encoding the stimulus category in a binary choice task,
#'    should be a factor with two levels, otherwise it will be transformed to
#'    a factor with a warning),
#' * \code{correct} (encoding whether the decision was correct; values in 0, 1),
#' * \code{sbj} (giving the subject ID; the models given in the second argument are fitted for each
#'   subject individually.
#' @param models `character` vector of models to be fit for each participant.
#' Implemented models: 'WEV', 'SDT', 'Noisy', 'PDA', and '2Chan'.
#' Alternatively, if `models="all"` (default), all implemented models will be fit.
#' @param var `character`. One of "constant" (default), "increasing", or "free" (will be
#' implemented), indicating how noise variances should be treated across conditions.
#' See Details for more information. Applies only to the models "SDT" and "WEV".
#' @param .parallel `logical`. Whether to parallelise the fitting over models and participant
#' (default: FALSE)
#' @param n.cores `integer`. Number of cores used for parallelization. If NULL (default), the available
#' number of cores -1 will be used.
#'
#' @return Gives data frame with rows for each model-participant combination and columns for the different parameters
#' as fitted result as well as additional information about the fit (`negLogLik` (for final parameters),
#' `k` (number of parameters), `N` (number of data rows), `BIC`, `AICc` and `AIC`)
#'
#'
#' @details The fitting involves a first grid search through an initial grid. Then the best 25 (\code{nAttempts})
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
#' @name fitConfModels
#' @import parallel
#' @importFrom stats dnorm pnorm optim integrate
# @importFrom pracma integral
#'
#' @references Rausch, M., Hellmann, S. & Zehetleitner, M. (2018). Confidence in
#'  masked orientation judgments is informed by both evidence and visibility.
#'  \emph{Atten Percept Psychophys} 80, 134–154. doi: 10.3758/s13414-017-1431-5
#'
#'
#' @examples
#' # 1. Generate data from two artificial participants
#' ###      d1,  d2, d3, cA3, cA2, cA1,theta, cB1, cB2, cB3)
#' p1 <- c(0.2, 0.5, 1.5,-1.0,-0.5,-0.2,  0.1, 0.3, 0.6, 0.9)
#' p2 <- c(0.1, 0.7, 1.9,-0.8,-0.4,-0.1,  0.0, 0.1, 0.4, 0.8)
#' D <- rbind(p1[1:3], p2[1:3])
#' thresholds <- rbind(c(-Inf, p1[4:10], Inf),
#'                     c(-Inf, p2[4:10], Inf))
#' df <- expand.grid(n = 1:30, stimulus = c(-1, 1), condition=c(1,2, 3))
#' data <- data.frame()
#' for (i in 1:2) { #for each participant
#'   x <- rnorm(nrow(df), mean=df$stimulus*D[i, df$condition])
#'   response1 <- as.numeric(cut(x,breaks=thresholds[i,], include.lowest = TRUE))
#'   data <- rbind(data, cbind(df, x, response1, participant=i))
#' }
#' data$response <- ifelse(data$response1<=4, -1, 1)
#' data$correct <- as.numeric(data$response==data$stimulus)
#' data$rating <- as.factor(ifelse(data$response==-1, 5-data$response1, data$response1-4))
#' data$stimulus <- as.factor(data$stimulus)
#' data$condition <- as.factor(data$condition)
#' table(data[data$correct==1,c("participant", "condition")])/60
#' head(data)
#'
#'
#' # 2. Use fitting function
#' \dontrun{
#'   # Fitting takes very long to run and uses multiple cores with this
#'   # call:
#'   fitConfModels(data, models=c("SDT", "WEV", "Noisy"),
#'                 .parallel=TRUE)
#' }
#'
#'

#' @export
fitConfModels <- function(data, models="all", #var="constant",
                             .parallel=FALSE, n.cores=NULL) {
  AllModels <- c('WEV', 'SDT','IG','ITGc', 'ITGcm', 'Noisy', 'PDA') # if you implement additional models, add them here!
  if (identical(models,"all")) models <- AllModels
  if (!all(models %in% AllModels)) {
    stop(paste(paste(setdiff(models, AllModels),collapse = " and "), " not implemented!"))
  }
  if (length(unique(models))<length(models)) {
    warning("Duplicate models are dropped")
    models <- unique(models)
  }
  if (is.null(data$condition)) data$condition <- factor(1)
  if (!is.factor(data$condition)) {
    data$condition <- factor(data$condition)
    warning("condition is transformed to a factor!")
  }
  if (!is.factor(data$stimulus)) {
    data$stimulus <- factor(data$stimulus)
    warning("stimulus is transformed to a factor!")
  }
  if (!is.factor(data$rating)) {
    data$rating <- factor(data$rating)
    warning("rating is transformed to a factor!")
  }
  nConds <- length(unique(data$condition))
  nRatings <- length(unique(data$rating))
  ## Define common names for the output to rbind all parameter fits together
  ## ToDo: Namen anpassen
  outnames <- c("model", "participant", "negLogLik", "N", "k", "BIC", "AICc", "AIC",
                paste("d", 1:nConds, sep=""),
                "theta", "w", "a", "sigma", "m_ratio",
                paste0("cA",1:(nRatings-1)),
                paste0("cB",1:(nRatings-1)))
  # This function will be called for every combination of participant and model
  call_fitfct <- function(X) {
    cur_model <- models[X[1]]
    cur_sbj <- X[2]
    participant <- NULL # to omit a note in R checks because of an unbound variable
    data_part <- subset(data, participant==cur_sbj)
    res <- fitConf(data_part, model = cur_model, var)
    res$model <- cur_model
    res$participant <- cur_sbj
    res[outnames[!(outnames %in% names(res))]] <- NA
    res <- res[,outnames]
    return(res)
  }

  # generate a list of fitting jobs to do and setup parallelization
  subjects <- unique(data$participant)
  nJobs <- length(models)*length(subjects)
  jobs <- expand.grid(model=1:length(models), sbj=subjects)
  if (.parallel) {
    listjobs <- list()
    for (i in 1:nrow(jobs)) {
      listjobs[[i]] <- c(model = jobs[["model"]][i], sbj = jobs[["sbj"]][i])
    }
    if (is.null(n.cores)) n.cores <- parallel::detectCores()-1

    cl <- makeCluster(type="SOCK", n.cores)
    clusterExport(cl, c("data",  "models","outnames", "call_fitfct", "var"),
                  envir = environment())
    # Following line ensures that the cluster is stopped even in cases of user
    # interrupt or errors
    on.exit(try(stopCluster(cl), silent = TRUE))
    res <- clusterApplyLB(cl, listjobs, fun=call_fitfct)
    stopCluster(cl)
  } else {
    res <- apply(X=jobs, 1, FUN=call_fitfct)
  }
  # bind list-outout together into data.frame
  res <- do.call(rbind, res)

  # finally, drop columns with unnecessary parameters
  res <- res[,apply(res, 2, function(X) any(!is.na(X)))]
  return(res)
}

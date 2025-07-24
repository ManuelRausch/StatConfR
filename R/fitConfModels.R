#' @title Fit several static confidence models to multiple participants

#' @description The `fitConfModels` function fits the parameters of several computational models of decision
#' confidence, in binary choice tasks,  specified in the `model` argument, to
#' different subsets of one data frame, indicated by different values in the column
#' `participant` of the `data` argument.
#' `fitConfModels` is a wrapper of the function \code{\link{fitConf}} and calls
#'  \code{\link{fitConf}} for every possible combination
#' of model in the `models` argument and sub-data frame of `data` for each value
#' in the `participant` column.
#' See Details for more information about the parameters.
#' Parameters are fitted using a maximum likelihood estimation method with a
#' initial grid search to find promising starting values for the optimization.
#' In addition, several measures of model fit (negative log-likelihood, BIC, AIC, and AICc)
#' are computed, which can be used for a quantitative model evaluation.

#' @param data  a `data.frame` where each row is one trial, containing following
#' variables:
#' * \code{diffCond} (optional; different levels of discriminability,
#'    should be a factor with levels ordered from hardest to easiest),
#' * \code{rating} (discrete confidence judgments, should be a factor with levels
#'    ordered from lowest confidence to highest confidence;
#'    otherwise will be transformed to factor with a warning),
#' * \code{stimulus} (stimulus category in a binary choice task,
#'    should be a factor with two levels, otherwise it will be transformed to
#'    a factor with a warning),
#' * \code{correct} (encoding whether the response was correct; should  be 0 for
#'    incorrect responses and 1 for correct responses)
#' * \code{participant} (some group ID, most often a participant identifier;
#'    the models given in the second argument are fitted to each subset of `data`
#'    determined by the different values of this column)
#' @param models `character`. The different computational models that should be
#'    fitted. Models implemented so far: 'WEV', 'SDT', 'GN', 'PDA', 'IG', 'ITGc',
#'    'ITGcm', 'logN', and 'logWEV'. Alternatively, if `model="all"` (default),
#'    all implemented models will be fit.
#' @param nInits `integer`. Number of initial values used for maximum likelihood optimization.
#' Defaults to 5.
#' @param nRestart `integer`. Number of times the optimization is restarted.
#' Defaults to 4.
#' @param .parallel `logical`. Whether to parallelize the fitting over models and participant
#' (default: FALSE)
#' @param n.cores `integer`. Number of cores used for parallelization. If NULL (default), the available
#' number of cores -1 will be used.

#' @return Gives `data.frame` with one row for each combination of model and
#' participant. There are different columns for the model, the participant ID, and one
#' one column for each estimated model parameter (parameters
#' not present in a specific model are filled with NAs).
#' Additional information  about the fit is provided in additional columns:
#' - `negLogLik` (negative log-likelihood of the best-fitting set of parameters),
#' - `k` (number of parameters),
#' - `N` (number of trials),
#' - `AIC` (Akaike Information Criterion; Akaike, 1974),
#' - `BIC` (Bayes information criterion; Schwarz, 1978),
#' - `AICc` (AIC corrected for small samples; Burnham & Anderson, 2002)
#' If length(models) > 1 or models == "all", there will be three additional columns:

#' @details
#' The provided `data` argument is split into subsets according to the values of
#' the `participant` column. Then for each subset and each model in the `models`
#' argument, the parameters of the respective model are fitted to the data subset.
#'
#' The fitting routine first performs a coarse grid search to find promising
#' starting values for the maximum likelihood optimization procedure. Then the best \code{nInits}
#' parameter sets found by the grid search are used as the initial values for separate
#' runs of the Nelder-Mead algorithm implemented in \code{\link[stats]{optim}}.
#' Each run is restarted \code{nRestart} times.
#'
#' ## Mathematical description of models
#'
#' The computational models are all based on signal detection theory (Green & Swets, 1966). It is assumed
#' that participants select a binary discrimination response \eqn{R} about a stimulus \eqn{S}.
#' Both \eqn{S} and \eqn{R} can be either -1 or 1.
#' \eqn{R} is considered correct if \eqn{S=R}.
#' In addition, we assume that there are \eqn{K} different levels of stimulus discriminability
#' in the experiment, i.e. a physical variable that makes the discrimination task easier or harder.
#' For each level of discriminability, the function fits a different discrimination
#' sensitivity parameter \eqn{d_k}. If there is more than one sensitivity parameter,
#' we assume that the sensitivity parameters are ordered such as \eqn{0 < d_1 < d_2 < ... < d_K}.
#' The models assume that the stimulus generates normally distributed sensory evidence \eqn{x} with mean \eqn{S\times d_k/2}
#' and variance of 1. The sensory evidence \eqn{x} is compared to a decision
#'  criterion \eqn{c} to generate a discrimination response
#' \eqn{R}, which is 1, if \eqn{x} exceeds \eqn{c} and -1 else.
#' To generate confidence, it is assumed that the confidence variable \eqn{y} is compared to another
#' set of criteria \eqn{\theta_{R,i}, i=1,2,...,L-1}, depending on the
#' discrimination response \eqn{R} to produce a \eqn{L}-step discrete confidence response.
#' The number of thresholds will be inferred from the number of steps in the
#' `rating` column of `data`.
#'  Thus, the parameters shared between all models are:
#' - sensitivity parameters \eqn{d_1},...,\eqn{d_K} (\eqn{K}: number of difficulty levels)
#' - decision criterion \eqn{c}
#' - confidence criterion \eqn{\theta_{-1,1}},\eqn{\theta_{-1,2}},
#' ..., \eqn{\theta_{-1,L-1}}, \eqn{\theta_{1,1}},  \eqn{\theta_{1,2}},...,
#' \eqn{\theta_{1,L-1}} (\eqn{L}: number of confidence categories available for confidence ratings)
#'
#' How the confidence variable \eqn{y} is computed varies across the different models.
#' The following models have been implemented so far:
#'
#' ### \strong{Signal detection rating model (SDT)}
#' According to SDT, the same sample of sensory
#' evidence is used to generate response and confidence, i.e.,
#' \eqn{y=x} and the confidence criteria span from the left and
#' right side of the decision criterion \eqn{c}(Green & Swets, 1966).
#'
#' ### \strong{Gaussian noise model (GN)}
#' According to the model, \eqn{y} is subject to
#' additive noise and assumed to be normally distributed around the decision
#' evidence value \eqn{x} with a standard deviation \eqn{\sigma}(Maniscalco & Lau, 2016).
#' \eqn{\sigma} is an additional free parameter.
#'
#' ### \strong{Weighted evidence and visibility model (WEV)}
#' WEV assumes that the observer combines evidence about decision-relevant features
#' of the stimulus with the strength of evidence about choice-irrelevant features
#' to generate confidence (Rausch et al., 2018). Thus, the WEV model assumes that \eqn{y} is normally
#' distributed with a mean of \eqn{(1-w)\times x+w \times d_k\times R} and standard deviation \eqn{\sigma}.
#' The standard deviation quantifies the amount of unsystematic variability
#' contributing to confidence judgments but not to the discrimination judgments.
#' The parameter \eqn{w} represents the weight that is put on the choice-irrelevant
#' features in the confidence judgment. \eqn{w} and \eqn{\sigma} are fitted in
#' addition to the set of shared parameters.
#'
#' ### \strong{Post-decisional accumulation model (PDA)}
#' PDA represents the idea of on-going information accumulation after the
#' discrimination choice (Rausch et al., 2018). The parameter \eqn{a} indicates the amount of additional
#' accumulation. The confidence variable is normally distributed with mean
#' \eqn{x+S\times d_k\times a} and variance \eqn{a}.
#' For this model the parameter \eqn{a} is fitted in addition to the shared
#' parameters.
#'
#' ### \strong{Independent Gaussian model (IG)}
#' According to IG, \eqn{y} is sampled independently
#' from \eqn{x} (Rausch & Zehetleitner, 2017). \eqn{y} is normally distributed with a mean of \eqn{a\times d_k} and variance
#' of 1 (again as it would scale with \eqn{m}). The additional parameter \eqn{m}
#' represents the amount of information available for confidence judgment
#' relative to amount of evidence available for the discrimination decision and can
#'  be smaller as well as greater than 1.
#'
#' ### \strong{Independent truncated Gaussian model: HMetad-Version (ITGc)}
#' According to the version of ITG consistent
#' with the HMetad-method (Fleming, 2017; see Rausch et al., 2023), \eqn{y} is sampled independently
#' from \eqn{x} from a truncated Gaussian distribution with a location parameter
#' of \eqn{S\times d_k \times m/2} and a scale parameter of 1. The Gaussian distribution of \eqn{y}
#' is truncated in a way that it is impossible to sample evidence that contradicts
#' the original decision: If \eqn{R = -1}, the distribution is truncated to the
#' right of \eqn{c}. If \eqn{R = 1}, the distribution is truncated to the left
#' of \eqn{c}. The additional parameter \eqn{m} represents metacognitive efficiency,
#' i.e., the amount of information available for confidence judgments relative to
#' amount of evidence available for discrimination decisions and  can be smaller
#' as well as greater than 1.
#'
#' ### \strong{Independent truncated Gaussian model: Meta-d'-Version (ITGcm)}
#' According to the version of the ITG consistent
#' with the original meta-d' method (Maniscalco & Lau, 2012, 2014; see Rausch et al., 2023),
#' \eqn{y} is sampled independently from \eqn{x} from a truncated Gaussian distribution with a location parameter
#' of \eqn{S\times d_k \times m/2} and a scale parameter
#' of 1. If \eqn{R = -1}, the distribution is truncated to the right of \eqn{m\times c}.
#' If \eqn{R = 1}, the distribution is truncated to the left of  \eqn{m\times c}.
#' The additional parameter \eqn{m} represents metacognitive efficiency, i.e.,
#' the amount of information available for confidence judgments relative to
#' amount of evidence available for the discrimination decision and  can be smaller
#' as well as greater than 1.
#'
#' ### \strong{Logistic noise model (logN)}
#' According to logN, the same sample
#' of sensory evidence is used to generate response and confidence, i.e.,
#' \eqn{y=x} just as in SDT (Shekhar & Rahnev, 2021). However, according to logN, the confidence criteria
#' are not assumed to be constant, but instead they are affected by noise drawn from
#' a lognormal distribution. In each trial, \eqn{\theta_{-1,i}} is given
#' by \eqn{c -  \epsilon_i}. Likewise,  \eqn{\theta_{1,i}} is given by
#' \eqn{c + \epsilon_i}. \eqn{\epsilon_i} is drawn from a lognormal distribution with
#' the location parameter
#' \eqn{\mu_{R,i}=log(|\overline{\theta}_{R,i}- c|) - 0.5 \times \sigma^{2}} and
#' scale parameter \eqn{\sigma}. \eqn{\sigma} is a free parameter designed to
#' quantify metacognitive ability. It is assumed that the criterion noise is perfectly
#' correlated across confidence criteria, ensuring that the confidence criteria
#' are always perfectly ordered. Because \eqn{\theta_{-1,1}}, ..., \eqn{\theta_{-1,L-1}},
#' \eqn{\theta_{1,1}}, ..., \eqn{\theta_{1,L-1}} change from trial to trial, they are not estimated
#' as free parameters. Instead, we estimate the means of the confidence criteria, i.e., \eqn{\overline{\theta}_{-1,1}, ...,
#' \overline{\theta}_{-1,L-1}, \overline{\theta}_{1,1}, ...  \overline{\theta}_{1,L-1}},
#' as free parameters.
#'
#' ### \strong{Logistic weighted evidence and visibility model (logWEV)}
#' logWEV is a combination of logN and WEV proposed by Shekhar and Rahnev (2023).
#' Conceptually, logWEV assumes that the observer combines evidence about decision-relevant features
#' of the stimulus with the strength of evidence about choice-irrelevant features (Rausch et al., 2018).
#' The model also assumes that noise affecting the confidence decision variable is lognormal
#'  in accordance with Shekhar and Rahnev (2021).
#' According to logWEV, the confidence decision variable is \eqn{y} is equal to
#' \eqn{y^*\times R}. \eqn{y^*} is sampled from a lognormal distribution with a location parameter
#'  of \eqn{(1-w)\times x\times R + w \times d_k} and a scale parameter of \eqn{\sigma}.
#'  The parameter \eqn{\sigma} quantifies the amount of unsystematic variability
#' contributing to confidence judgments but not to the discrimination judgments.
#' The parameter \eqn{w} represents the weight that is put on the choice-irrelevant
#' features in the confidence judgment. \eqn{w} and \eqn{\sigma} are fitted in
#' addition to the set of shared parameters.

#' @author
#' Sebastian Hellmann, \email{sebastian.hellmann@tum.de}\cr
#' Manuel Rausch, \email{manuel.rausch@hochschule-rhein-waal.de}

# unlike for the other tags, the references are formatted more nicely if each reference is tagged seperately
#' @references Akaike, H. (1974). A New Look at the Statistical Model Identification. IEEE Transactions on Automatic Control, AC-19(6), 716–723.doi: 10.1007/978-1-4612-1694-0_16\cr
#' @references Burnham, K. P., & Anderson, D. R. (2002). Model selection and multimodel inference: A practical information-theoretic approach. Springer.\cr
#' @references Fleming, S. M. (2017). HMeta-d: Hierarchical Bayesian estimation of metacognitive efficiency from confidence ratings. Neuroscience of Consciousness, 1, 1–14. doi: 10.1093/nc/nix007\cr
#' @references Green, D. M., & Swets, J. A. (1966). Signal detection theory and psychophysics. Wiley.\cr
#' @references Maniscalco, B., & Lau, H. (2012). A signal detection theoretic method for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430.\cr
#' @references Maniscalco, B., & Lau, H. C. (2014). Signal Detection Theory Analysis of Type 1 and Type 2 Data: Meta-d’, Response- Specific Meta-d’, and the Unequal Variance SDT Model. In S. M. Fleming & C. D. Frith (Eds.), The Cognitive Neuroscience of Metacognition (pp. 25–66). Springer. doi: 10.1007/978-3-642-45190-4_3\cr
#' @references Maniscalco, B., & Lau, H. (2016). The signal processing architecture underlying subjective reports of sensory awareness. Neuroscience of Consciousness, 1, 1–17. doi: 10.1093/nc/niw002\cr
#' @references Rausch, M., Hellmann, S., & Zehetleitner, M. (2018). Confidence in masked orientation judgments is informed by both evidence and visibility. Attention, Perception, and Psychophysics, 80(1), 134–154. doi: 10.3758/s13414-017-1431-5\cr
#' @references Rausch, M., Hellmann, S., & Zehetleitner, M. (2023). Measures of metacognitive efficiency across cognitive models of decision confidence. Psychological Methods. doi: 10.31234/osf.io/kdz34\cr
#' @references Rausch, M., & Zehetleitner, M. (2017). Should metacognition be measured by logistic regression? Consciousness and Cognition, 49, 291–312. doi: 10.1016/j.concog.2017.02.007\cr
#' @references Schwarz, G. (1978). Estimating the dimension of a model. The Annals of Statistics, 6(2), 461–464. doi: 10.1214/aos/1176344136\cr
#' @references Shekhar, M., & Rahnev, D. (2021). The Nature of Metacognitive Inefficiency in Perceptual Decision Making. Psychological Review, 128(1), 45–70. doi: 10.1037/rev0000249\cr
#' @references Shekhar, M., & Rahnev, D. (2023). How Do Humans Give Confidence? A Comprehensive Comparison of Process Models of Perceptual Metacognition. Journal of Experimental Psychology: General. doi:10.1037/xge0001524\cr

#' @examples
#' # 1. Select two subjects from the masked orientation discrimination experiment
#' data <- subset(MaskOri, participant %in% c(1:2))
#' head(data)
#'
#' # 2. Fit some models to each subject of the masked orientation discrimination experiment
#' \donttest{
#'   # Fitting several models to several subjects takes quite some time
#'   # (about 10 minutes per model fit per participant on a 2.8GHz processor
#'   # with the default values of nInits and nRestart).
#'   # If you want to fit more than just two subjects,
#'   # we strongly recommend setting .parallel=TRUE
#'   Fits <- fitConfModels(data, models = c("SDT", "ITGc"), .parallel = FALSE)
#' }

#' @import parallel
#' @importFrom stats dnorm pnorm qnorm optim integrate plnorm

#' @export
fitConfModels <- function(data, models="all",
                         # diffCond = NULL, stimulus = NULL, correct = NULL, rating = NULL,
                          nInits = 5, nRestart = 4,
                          .parallel=FALSE, n.cores=NULL) {
  AllModels <- c('WEV', 'SDT', 'IG', 'ITGc',
                 'ITGcm', 'GN', 'PDA', 'logN', 'logWEV', 'RCE') # if you implement additional models, add them here!
  if (identical(models,"all")) models <- AllModels
  if (!all(models %in% AllModels)) {
    stop(paste(paste(setdiff(models, AllModels),collapse = " and "), " not implemented!"))
  }
  if (length(unique(models))<length(models)) {
    warning("Duplicate models are dropped")
    models <- unique(models)
  }
  # if (!is.null(diffCond)) data$diffCond <- data[,diffCond]
  # if (!is.null(stimulus)) data$stimulus <- data[,stimulus]
  # if (!is.null(correct)) data$correct <- data[,correct]
  # if (!is.null(rating)) data$rating <- data[,rating]
  if (is.null(data$diffCond)) data$diffCond <- factor(1)
  if (!is.factor(data$diffCond)) {
    data$diffCond <- factor(data$diffCond)
    warning("diffCond is transformed to a factor!")
  }
  if(length(unique(data$stimulus)) != 2) {
    stop("There must be exacltly two different possible values of stimulus")
  }

  if (!is.factor(data$stimulus)) {
    data$stimulus <- factor(data$stimulus)
    warning("stimulus is transformed to a factor!")
  }
  if (!is.factor(data$rating)) {
    data$rating <- factor(data$rating)
    warning("rating is transformed to a factor!")
  }
  if(!all(data$correct %in% c(0,1))) stop("correct should be 1 or 0")

  nConds <- length(unique(data$diffCond))
  nRatings <- length(unique(data$rating))
  ## Define common names for the output to rbind all parameter fits together
  ## ToDo: Namen anpassen
  outnames <- c("model", "participant", "negLogLik", "N", "k", "BIC", "AICc", "AIC",
                paste("d_", 1:nConds, sep=""),"c",
                paste0("theta_minus.",(nRatings-1):1),
                paste0("theta_plus.",1:(nRatings-1)),
                paste0("M_theta_minus.",(nRatings-1):1),
                paste0("M_theta_plus.",1:(nRatings-1)),
                "b", "m", "sigma", "w"
  )
  # This function will be called for every combination of participant and model
  call_fitfct <- function(X) {
    cur_model <- models[X[1]]
    cur_sbj <- X[2]
    participant <- NULL # to omit a note in R checks because of an unbound variable
    data_part <- subset(data, participant==cur_sbj)
    res <- fitConf(data_part, model = cur_model,  nInits = nInits, nRestart = nRestart)
    res$model <- cur_model
    res$participant <- cur_sbj
    res[outnames[!(outnames %in% names(res))]] <- NA
    res <- res[,outnames]
    return(res)
  }

  # generate a list of fitting jobs to do and setup parallelization
  no_sbj_column <- FALSE
  if (is.null(data$participant)) {
    data$participant <- 1
    no_sbj_column <- TRUE
  }
  subjects <- unique(data$participant)
  nJobs <- length(models)*length(subjects)
  jobs <- expand.grid(model=1:length(models), sbj=subjects)
  if (.parallel) {
    listjobs <- list()
    for (i in 1:nrow(jobs)) {
      listjobs[[i]] <- c(model = jobs[["model"]][i], sbj = jobs[["sbj"]][i])
    }
    if (is.null(n.cores)) n.cores <- min(nJobs, detectCores()-1)

    cl <- makeCluster(type="SOCK", n.cores)
    clusterExport(cl, c("data",  "models", "outnames", "call_fitfct", "nInits", "nRestart"),
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
  if (no_sbj_column) res$participant <- NULL
  return(res)
}

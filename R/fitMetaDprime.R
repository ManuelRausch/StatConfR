#' Fits meta-d′ and meta-d′/d′ ratios for data from one or several subjects
#'
#' This function computes meta-d′ and meta-d′/d′ for each  participant in the \code{data}, respectively.
#' @param data  a `data.frame` where each row is one trial, containing following
#' variables:
#' * \code{rating} (discrete confidence judgments, should be given as factor;
#'    otherwise will be transformed to factor with a warning),
#' * \code{stimulus} (stimulus category in a binary choice task,
#'    should be a factor with two levels, otherwise it will be transformed to
#'    a factor with a warning),
#' * \code{correct} (encoding whether the response was correct; should  be 0 for incorrect responses and 1 for correct responses)
#' * \code{participant} (giving the subject ID; the models given in the second argument are fitted for each
#'   subject individually.
#' @param model `character` of length 1. Either "ML" to use the original model
#' specification by Maniscalco and Lau (2012,  2014) or "F" to use the model
#' specification by Fleming (2017)'s HmetaD method.  Defaults to "ML"
#' @param nInits `integer`. Number of initial values used for maximum likelihood optimization.
#' Defaults to 5.
#' @param nRestart `integer`. Number of times the optimization is restarted.
#' Defaults to 3.
#' @param .parallel `logical`. Whether to parallelize the fitting over models and participant
#' (default: FALSE)
#' @param n.cores `integer`. Number of cores used for parallelization. If NULL (default), the available
#' number of cores -1 will be used.
#'
#' @return Gives data frame with rows for each participant and columns dprime, c, metaD, and Ratio
#' - dprime is the discrimination sensitivity index d′, calculated using a standard SDT formula
#' - c is the discrimination bias c, calculated using a standard SDT formula
#' - metaD is meta-d′, discrimination sensitivity estimated from confidence judgments conditioned on the response
#' - Ratio is meta-d′/d′, a quantity usually referred to as metacognitive efficiency.
#'
#' @details
#' The function computes meta-d′ and meta-d′/d′ either using the
#' hypothetical signal detection model assumed by Maniscalco and Lau (2012, 2014)
#' or the one assumed by Fleming (2014). The fitting routine first performs a coarse grid search to find promising
#' starting values for the maximum likelihood optimization procedure. Then the best \code{nInits}
#' parameter sets found by the grid search are used as the initial values for separate
#' runs of the Nelder-Mead algorithm implemented in \code{\link[stats]{optim}}.
#' Each run is restarted \code{nRestart} times. Warning: meta-d′/d′
#' is only guaranteed to be unbiased from discrimination sensitivity, discrimination
#' bias, and confidence criteria if the data is generated according to the
#' independent truncated Gaussian model (see Rausch et al., 2023).
#'
#' @md
#'
#' @author Manuel Rausch, \email{manuel.rausch@@hochschule-rhein-waal.de}
#'
#' @name fitMetaDprime
#'
#' @import parallel
#' @importFrom stats dnorm pnorm pnorm optim integrate
#'
#' @references Fleming, S. M. (2017). HMeta-d: Hierarchical Bayesian estimation of metacognitive efficiency from confidence ratings. Neuroscience of Consciousness, 1, 1–14. doi: 10.1093/nc/nix007
#' @references Maniscalco, B., & Lau, H. (2012). A signal detection theoretic method for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430.
#' @references Maniscalco, B., & Lau, H. C. (2014). Signal Detection Theory Analysis of Type 1 and Type 2 Data: Meta-d’, Response- Specific Meta-d’, and the Unequal Variance SDT Model. In S. M. Fleming & C. D. Frith (Eds.), The Cognitive Neuroscience of Metacognition (pp. 25–66). Springer. doi: 10.1007/978-3-642-45190-4_3
#' @references Rausch, M., Hellmann, S., & Zehetleitner, M. (2023). Measures of metacognitive efficiency across cognitive models of decision confidence (Preprint). PsyArXiv. doi: 10.31234/osf.io/kdz34
#'
#' @examples
#' # 1. Select two subject from the masked orientation discrimination experiment
#' data <- subset(MaskOri, participant %in% c(1:2))
#' head(data)
#'
#' # 2. Fit meta-d′/d′ for each subject in data
#' MetaDs <- fitMetaDprime(data, model="F", .parallel = FALSE)
#'


#' @export
fitMetaDprime <- function(data, model="ML",  nInits = 5, nRestart = 3,
                          .parallel=FALSE, n.cores=NULL) {
  if (! model %in% c("ML", "F")) {
    stop("model must be either 'ML' or 'F'")
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
  if(!any(data$correct == 0)) stop("There should be at least one erroneous response")
  if(!any(data$correct == 1)) stop("There should be at least one correct response")

  nRatings <- length(unique(data$rating))
  ## Define common names for the output to rbind all parameter fits together
  ## ToDo: Namen anpassen
  outnames <- c("model", "participant", "dprime", "c", "metaD", "Ratio")
  # This function will be called for every combination of participant and model
  call_fitfct <- function(X) {
    cur_model <- model[X[1]]
    cur_sbj <- X[2]
    participant <- NULL # to omit a note in R checks because of an unbound variable
    data_part <- subset(data, participant==cur_sbj)
    res <- int_fitMetaDprime(ratings=data_part$rating,
                             stimulus=data_part$stimulus, correct = data_part$correct,
                             ModelVersion = cur_model,  nInits = nInits, nRestart = nRestart)
    res$model <- cur_model
    res$participant <- cur_sbj
    res[outnames[!(outnames %in% names(res))]] <- NA
    res <- res[,outnames]
    return(res)
  }

  # generate a list of fitting jobs to do and setup parallelization
  subjects <- unique(data$participant)
  nJobs <- length(model)*length(subjects)
  jobs <- expand.grid(model=1:length(model), sbj=subjects)
  if (.parallel) {
    listjobs <- list()
    for (i in 1:nrow(jobs)) {
      listjobs[[i]] <- c(model = jobs[["model"]][i], sbj = jobs[["sbj"]][i])
    }
    if (is.null(n.cores)) n.cores <- detectCores() - 1

    cl <- makeCluster(type="SOCK", n.cores)
    clusterExport(cl, c("data",  "model","outnames", "call_fitfct", "nInits", "nRestart"),
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

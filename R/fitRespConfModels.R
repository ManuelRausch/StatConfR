#' fit several static confidence models to multiple participants
#'
#' @param data data.frame containing the data to fit the models to. Should contain
#' following columns: participant, stimulus, correct, and rating
#' @param models character. Vector of models to be fit for each participant.
#' Implemented models: 'WEV', 'SDT', 'Noisy', 'PDA', and '2Chan'.
#' @param .parallel logical. Whether to parallelise the fitting over models and participant
#' @param n.cores integer. Number of cores used for parallelisation. If NULL, the available
#' number of cores -1 will be used.
#'
#' @return data.frame with one row per combination of participant and model and
#' columns for participant, model, the fitted parameters and additional information
#' about the fit (BIC, AIC,...)
#' @export
#'
#' @examples
#' #ToDo
#'
#' @import parallel
fitMultiRespConf <- function(data, models,
                             .parallel=FALSE, n.cores=NULL) {
  if (!all(models %in% c('WEV', 'SDT','2Chan','Noisy', 'PDA'))) {
    stop("At least one model in models argument is not implemented!")
  }
  if (length(unique(models))<length(models)) {
    warning("There were duplicate models, which are dropped")
    models <- unique(models)
  }
  ## Define common names for the output to rbind all parameter fits together
  ## ToDo: Namen anpassen
  outnames <- c("model", "sbj", "negLogLik", "N", "k", "BIC", "AICc", "AIC", "fixed",
                "t0", "st0",
                paste("v", 1:nConds, sep=""),
                paste("thetaLower", 1:(nRatings-1), sep=""),
                paste("thetaUpper", 1:(nRatings-1), sep=""),
                "wrt", "wint", "wx", "b", "a",
                "z", "sz", "sv", "tau", "w", "svis", "sigvis", "lambda")
  # This function will be called for every participant and model combination
  call_fitfct <- function(X) {
    cur_model <- models[X[1]]
    cur_sbj <- X[2]
    participant <- NULL # to omit a note in R checks because of an unbound variable
    data_part <- subset(data, participant==cur_sbj)
    res <- fitRespConf(data_part, model = cur_model)
    res$model <- cur_model
    res$sbj <- cur_sbj
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
    clusterExport(cl, c("data",  "models","outnames", "call_fitfct"), envir = environment())
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

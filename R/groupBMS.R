#' @title Group-level Bayesian Model Comparison

#' @description \code{groupBMS} performs a Bayesian model comparison on the group level based on marginal
#' likelihoods (alias model evidence) obtained from model fits on the level of single subjects.
#' It compares three models of the distribution of model probabilities (see Daunizeau et al., 2014; Rigoux at al., 2014; Stephan et al., 2009):
#' * The random effects model treats the generative model as random effect that could differ between subjects. Model probabilities are assumed to follow a Dirichelet distribution in the population.
#' * The fixed effects model, according to which one single model generated the data of all subjects, but subjects can still different from each other due to different model parameters.
#' * The null model, according to which  model frequencies are a priori assumed to be fixed to 1/K, where K is the number of models to be compared.
#' \code{groupBMS} can be used with the output of \code{\link{fitConfModels}},
#' i.e. a data frame with information criteria for different models and subjects,
#' using one of the information criteria to approximate the model evidence.

#' @param fits a data frame as returned by \code{\link{fitConfModels}}.
#' Should contain a column `model`indicating the model name, a column
#' `participant`, indicating the grouping
#' structure of the data (e.g. participant IDs), and a column with the name
#' given by the `measure` argument containing the values of
#' the information criterion that should be used to approximate model evidence.
#'
#' @param measure the name of the the information criterion
#' to approximate model evidence. The available options are 'BIC', 'AIC', and 'AICc'.
#'
#' @param opts a list with options for the Variational Bayes algorithm to estimate
#' the parameter of the Dirichlet distribution. Following values can  be provided:
#' * \code{maxiter} the maximum number of iterations (Default: 200)
#' * \code{tol} the tolerance for changes in the free energy approximation
#'                 to stop the algorithm, if abs(FE(i+1)-FE(i))<tol the algorithm
#'                 is stopped (Default: 1e-4)
#' * \code{eps} The number to substitute values of 0 in calls to log (Default: 1e-32)
#'
#' @return \code{groupBMS} returns a list with two entries:
#' * `model_weights`: a dataframe with rows for each model and following columns:
#'    + `alpha`, the alpha parameter of the Dirichlet posterior over model probabilities in the population according to the random effects model,
#'    + `r`, the mean probabilities of each model in the population according to the random effects model,
#'    + `ep`, the exceedance probability for each model,
#'    + `pep`, the protected exceedance probabilities for each model, and
#'    + `fx_prop`, the posterior model probabilities if a fixed true model is assumed in the population).
#' * `summary_stats`: a vector giving statistics for the Bayesian model comparison analyses:
#'    + `bor_re_null`, Bayesian omnibus risk for the comparison of the random effects model against the null model,
#'    + `bor_fixed_null`, Bayesian omnibus risk for the comparison of the fixed effects model against the null model,
#'    + `bor_re_fixed`, Bayesian omnibus risk for the comparison of the random effects model against the fixed effects model
#'    + `FE_re`, estimated Free Energy according to the random effects model,
#'    + `FE_fixed` Free Energy according to the fixed effects model,
#'    + `FE_null` Free Energy according to the null model.
#'
#' @details The function `groupBMS` can be used for Bayesian model selection on
#' the group level following the approach of Rigoux et al. (2014).
#' It therefore aggregates the available subject-level model evidences
#' according to a random effects model, which assumes that model probabilities
#' are distributed in the population according to a Dirichlet distribution.
#' The parameters of the Dirichlet distribution are determined using a variational Bayes algorithm.
#' `groupBMS` provides three statistics associated with the random effects model:
#' * estimated model frequencies in the population,
#' * The exceedance probability  represents the probability that a specific
#' model is more probable than any other model (Stephan et al., 2009).
#' * The Protected exceedance probability represents the probability that a specific model
#'  is more probable than any other model for a random subject out of the population
#'  scaled by the Bayesian omnibus risk of assuming a random effects model over
#'  a null effect model, in which all models have the same probability (Rigoux et al., 2014).
#' The function is an R implementation of the VBA_groupBMC function of the VBA toolbox for MATLAB (Daunizeau, et al., 2014).
#'
#' @references Rigoux, L., Stephan, K. E., Friston, K. J., & Daunizeau, J. (2014). Bayesian model selection for group studies - revisited. \emph{NeuroImage}, 84, 971–985. doi: 10.1016/j.neuroimage.2013.08.065
#' @references Daunizeau, J., Adam, V., & Rigoux, L. (2014). Vba: A probabilistic treatment of nonlinear models for neurobiological and behavioural data. \emph{PLOS Computational Biology}, 10(1), e1003441. doi: 10.1371/journal.pcbi.1003441
#' @references Stephan, K. E., Penny, W. D., Daunizeau, J., Moran, R. J., & Friston, K. J. (2009). Bayesian model selection for group studies. Neuroimage, 46(4), 1004-1017.

#' @author
#' Sebastian Hellmann, \email{sebastian.hellmann@tum.de}\cr
#' Manuel Rausch, \email{manuel.rausch@ku.de}
#'
#' @examples
#' # Define a data frame with information criteria from model fits
#' # this is a sub-data.frame from an output of fitConfModels with
#' # 5 subjects, two models and rounded information criteria)
#' fitted_par <- data.frame(
#'   participant = rep(1:5, each=2),
#'   model = rep(c("WEV", "ITGcm"), 5),
#'   BIC = c(5360, 5550, 3773, 3963, 3441, 3503, 2613, 2706, 3566, 3695),
#'   AIC = c(5274, 5469, 3686, 3883, 3355, 3422, 2527, 2625, 3479, 3614),
#'   AICc = c(5275, 5469, 3687, 3883, 3356, 3422, 2527, 2625, 3480, 3614))
#' # Conduct group-level Bayesian model selection based on BIC
#' ModelComp <- groupBMS(fitted_par, measure="BIC")
#' ModelComp
#'
#' @importFrom stats rgamma
#' @export
groupBMS <- function(fits, measure = "AICc", opts=list()) {
  if (!measure %in% c("BIC", "AIC", "AICc"))
    stop(paste0("measure=", measure, " should be BIC, AIC, or AICc"))
  models <- sort(unique(fits$model))
  sbj_col <- c("sbj", "subject", "participant")
  sbj_col <- sbj_col[which(sbj_col %in% names(fits))]
  participants <- unique(fits[[sbj_col]])
  n <- length(participants)
  K <- length(models)
  fits <- fits[order(fits[[sbj_col]], fits[["model"]]),]
  mlp <- matrix(NA, ncol=n, nrow=K)
  for (i in 1:n) mlp[,i] <- fits[(i-1)*K + 1:K, measure]
  mlp <- - mlp/(2)
  rownames(mlp) <- models
  out <- .groupBMS(mlp, opts)
  return(out)
}

.groupBMS <- function(mlp, opts=list()) {
  if (is.null(opts$maxiter)) opts$maxiter <- 200 # max iterations for alpha-optimization
  if (is.null(opts$tol)) opts$tol <- 1e-4  # convergence criterion for alpha
  if (is.null(opts$eps)) opts$eps <- 1e-32   # replace 0s with eps in log-calls

  K <- nrow(mlp)
  n <- ncol(mlp)
  alpha0 <- rep(1, K)

  modelprobs_Null <- function(x) { # function for each n (x is a K-vector)
    g <- x - max(x)
    g <- exp(g)/sum(exp(g))
    res <- sum(g*(x - log(g+opts$eps) - log(K)))
    return(res)
  }

  FreeEnergyNull <- sum(apply(mlp, 2, modelprobs_Null))

  # derive probabilities and free energy of the 'fixed-effect' analysis
  ss <-  rowSums(mlp) + log(1/K)
  logz <- ss - max(ss)
  z <-  exp(logz)/sum(exp(logz))
  fixed_effects_postprobs  <- z
  FreeEnergy_fixed <- sum(z * ss) - sum(z * log(z+opts$eps))

  alpha <- alpha0
  FreeEnergy <- 0
  #Full_FEs <- NULL
  #ln_p_y_mk # log(p(y_n | m_nk)) = - BIC/2
  for (i in 1:opts$maxiter) {
    digamma_stuff <- digamma(alpha) - digamma(sum(alpha))
    logu_nk <- mlp + matrix(digamma_stuff, nrow=K, ncol=n, byrow=FALSE)
    ColMin <- apply(logu_nk, 2, min)
    logu_nk <- sweep(logu_nk, 2, ColMin)
    u_nk <- exp(logu_nk)

    u_nk[is.infinite(u_nk)] <- 1e+64
    model_sums <- colSums(u_nk)
    beta_nk <- sweep(u_nk, 2, model_sums, FUN="/")# z_nk bei Rigoux et al. (2014)
    # posterior.r in VBA toolbox

    Old_FreeEnergy <- FreeEnergy
    Sqf <- sum(lgamma(alpha)) - lgamma(sum(alpha)) - sum((alpha-1)*digamma_stuff)
    Sqm <- - sum(beta_nk * log(beta_nk + opts$eps))
    ELJ <- lgamma(sum(alpha0)) - sum(lgamma(alpha0)) + sum((alpha0-1)*digamma_stuff)
    ELJ <- ELJ + sum( beta_nk* (logu_nk))
    FreeEnergy <- Sqf + Sqm + ELJ

    beta_k <- rowSums(beta_nk)
    alpha <- alpha0 + beta_k
    #Full_FEs[i] <- FreeEnergy

    if (abs(Old_FreeEnergy-FreeEnergy)< opts$tol) { # equivalent to matlab (VBA-Toolbox TolFun = F(it)-F(it-1))
      # MaxIter in VBA toolbox is 32!
      #warning("Finished because of small update")
      break
    }
  }
  if (i == opts$maxiter) warning("Maximum number of iterations reached; alpha may have not converged!")

  ### Simulate exceedance probability
  simulate_ep <- function(n=100,alpha){
    # For dirichlet-rv we would scale Gamma-RV to a max of 1,
    # but we are interested in the index of the max only, so this it not
    # necessary
    K <- length(alpha)
    if (!is.null(names(alpha))) {
      res <- sapply(1:n, function(x) names(alpha)[which.max(rgamma(K, shape=alpha))])
      res <- factor(res, levels=names(alpha))
      res <- as.vector(table(res))/n
      names(res) <- names(alpha)
    } else {
      res <- sapply(1:n, function(x) which.max(rgamma(K, shape=alpha)))
      res <- as.vector(table(res))/n
    }
    return(res)
  }
  ep <- simulate_ep(1e+4, alpha)
  bor = 1/(1+exp(FreeEnergy-FreeEnergyNull))
  pep <- ep * (1 - bor) + bor / K

  bor_fixed <-     1/(1+exp(FreeEnergy_fixed-FreeEnergyNull))
  bor_re_fixed <-  1/(1+exp(FreeEnergy-FreeEnergy_fixed))

  model_weights <- matrix(c(alpha, alpha/sum(alpha), ep, pep, fixed_effects_postprobs), nrow=K)
  rownames(model_weights) <- rownames(mlp)
  colnames(model_weights) <- c("alpha", "r", "ep", "pep", "fx_prop")

  out <- list()
  out$model_weights <- model_weights
  out$summary_stats <- c(bor_re_null = bor,
                         bor_fixed_null = bor_fixed,
                         bor_re_fixed = bor_re_fixed,
                         FE_re = FreeEnergy,
                         FE_null=FreeEnergyNull,
                         FE_fixed=FreeEnergy_fixed)
  return(out)
}

#' Bayesian Model Comparison

#' \code{compareBayes} performs a Bayesian model comparison based on marginal
#' likelihoods (alias model evidence), given for different models across different
#' subject on a group level using a fixed effects model and a random effects model
#' on the distribution of model probabilities (see Rigoux at al., 2014; Daunizeau et al., 2014)
#' `compareBayes` can be used with the output of \code{\link{fitConfModels}}, i.e. a data frame with information
#' criteria for different models and subjects, using a information criterion to
#' approximate the model evidence.
#' \code{summaryCompareBayes} p
#' 

#' @param fits a data frame as returned by \code{\link{fitRTConfModels}}.
#' Should contain a column `model`indicating the model name, a column
#' `subject` (alternatively `sbj` or `participant`) indicating the grouping
#' structure of the data, and a column with the name given by the `measure`
#' argument containing the values of the information criterion that should be
#' used to approximate model evidence.
#' @param measure the name of the column indicating the information criterion
#' to approximate model evidence. For outputs of \code{\link{fitRTConfModels}},
#' the available measures are 'BIC', 'AIC', and 'AICc'. Any other approximation
#' for the model evidence may be used, the measure is transferred to log model
#' evidence by taking -measure/2.
#' @param opts a list with options for the iteration algorithm to estimate
#' the parameter of the Dirichlet distribution. Following values may be provided:
#' * \code{maxiter} the maximum number of iterations (Default: 200)
#' * \code{tol} the tolerance for changes in the free energy approximation
#'                 to stop the algorithm, if abs(FE(i+1)-FE(i))<tol the algorithm
#'                 is stopped (Default: 1e-4)
#' * \code{eps} The number to substitute values of 0 in calls to log (Default: 1e-32)
#' 
#' #' @return a matrix with rows for each model (row names indicate the
#'                     model names for `group_BMS_fits` and for `group_BMS` if
#'                     row names are available in `mlp`), and following columns:
#'                     `alpha` (the alpha parameter of the Dirichlet posterior
#'                     over model probabilities in the population), `r` (the
#'                     mean probabilities of each model in the population), `ep`
#'                     and `pep` (exceedance and protected exceedance
#'                     probabilities for each model), and `fx_prop` (the
#'                     posterior model probabilities if a fixed true model is
#'                     assumed in the population).                
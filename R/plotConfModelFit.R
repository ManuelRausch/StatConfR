#' @title Plot the prediction of fitted parameters of one model of confidence over the corresponding real data
#'
#' @description The `plotConfModelFit` function plots the predicted distribution of discrimination responses
#' and confidence ratings created from a `data.frame` of parameters obtaind from \code{\link{fitConfModels}}
#' and overlays the predicted distribution over the data to which the model parameters were fitted.
#'
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
#'
#' @param fitted_pars a `data.frame` with one row for each participant and model parameters in different columns.
#' fitted_pars also may contain a column called  `model`specifying the model to be visualized.
#' If there is no model column in data or if there are multiple models in fitted_pars,
#' it is necessary to specify the model argument.
#'
#' @param model `character`. See \code{\link{fitConf}} for all available models
#'
#' @return a `data.frame` with one row for each combination of model and
#' participant. There are different columns for the model, the participant ID, and one
#' one column for each estimated model parameter (parameters
#' not present in a specific model are filled with NAs)
#'
#' #' @examples
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
#'   Fits <- fitConfModels(data, models = "ITGc", .parallel = FALSE)
#' }
#' # 3. Plot the predicted probabilies based on model and fitted parameter over the observed relative frequencies.
#'
#' \donttest{
#'   # Fitting several models to several subjects takes quite some time
#'   # (about 10 minutes per model fit per participant on a 2.8GHz processor
#'   # with the default values of nInits and nRestart).
#'   # If you want to fit more than just two subjects,
#'   # we strongly recommend setting .parallel=TRUE
#'   myPlottedFit <- plotConfModelFit(data, Fits)
#'   myPlottedFit
#' }
#' @import ggplot2
#' @importFrom plyr ddply transform summarise
#' @importFrom Rmisc summarySEwithin
#'
#' @export
plotConfModelFit <- function(data, fitted_pars, model = NULL){
  if(is.null(model)){
    if("model" %in% colnames(fitted_pars)){
      if(length(unique(fitted_pars$model))==1){
        model <- unique(fitted_pars$model)
      } else {
        stop("Please use the model argument to specify which model should be used")
      }
    }else {
      stop("Please specify which model should be used")
    }
  }
  if (is.null(data$diffCond)) data$diffCond <- factor(1)
  if (!is.factor(data$diffCond)) {
    data$diffCond <- factor(data$diffCond)
  }
  if(length(unique(data$stimulus)) != 2) {
    stop("There must be exacltly two different possible values of stimulus")
  }

  if (!is.factor(data$stimulus)) {
    data$stimulus <- factor(data$stimulus)
  }
  if (!is.factor(data$rating)) {
    data$rating <- factor(data$rating)
  }
  if(!all(data$correct %in% c(0,1))) stop("correct should be 1 or 0")

  myColor <- switch(model, 'GN' = 1, 'IG' = 2, 'ITGc'  = 3, 'ITGcm' = 4, 'logN' = 5,
                    'logWEV' = 6,'PDA' = 7,  'WEV' = 8, 'SDT' = 9) # models are color coded

  # 1. First aggregate on the level of subjects

  AggDist <-
    plyr::ddply(data,
          ~  diffCond * rating * stimulus * correct * participant, #,
          plyr::summarise, p = length(rating),  .drop=FALSE)

  AggDist <- plyr::ddply(AggDist, ~ diffCond * stimulus,
                         transform, N = sum(p))
  AggDist$p <- AggDist$p / AggDist$N


  AggDist$rating <- as.numeric(AggDist$rating)
  AggDist$correct <-
    factor(AggDist$correct)
  levels(AggDist$correct)  <-
    c("incorrect", "correct")

   # 2. aggregate across subjects
  AggDist <-
    Rmisc::summarySEwithin(AggDist , measurevar = "p",
                           withinvars = c("diffCond", "correct", "rating", "stimulus"),
                           idvar = "participant",
                           na.rm = TRUE, .drop = TRUE)
  AggDist$rating <- as.numeric(AggDist$rating)
  levels(AggDist$stimulus) <- c("S = -1", "S = 1")
  levels(AggDist$diffCond) <- paste("K =", as.numeric(levels(AggDist$diffCond)))

  # 3) create a plot with the observed data
  PlotObsVsPred  <-
    ggplot(AggDist,
           aes(x=rating, y = p)) +
    facet_grid(diffCond ~ stimulus+correct) +
    geom_bar(stat="identity",
             fill ="white", color = "black") +
    geom_errorbar(aes(ymin=p-se, ymax=p+se), width=0) +
    xlab("Confidence rating") +
    ylab("probability") +
    theme(strip.text.y = element_text(angle=0)) +
    theme_minimal()


  # create a plot from the observed data


  PlotObsVsPred
}


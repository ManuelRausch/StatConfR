#' Get Normal Noise Information
#' 
#' Determine the transmitted information of a normal noise classifier. Either
#' accuracies or sensitivities have to be specified to determine the
#' corresponding transmitted information. If both are specified (not advised
#' due to ambiguity), accuracies are prioritizied.
#' 
#' @param accuracies Vector of accuracies
#' @param sensitivities Vector of sensitivities
#' 
#' @return Transmitted information of a normal noise classifier
get_normal_noise_information <- function(accuracies    = NULL,
                                         sensitivities = NULL)
{
  if (is.null(sensitivities) & is.null(accuracies))
    stop("Specify accuracies or sensitivities")

  if (is.null(sensitivities))
    sensitivities <- transform_normal_accuracy_to_sensitivity(accuracies)

  if (is.null(accuracies))
    accuracies <- transform_normal_sensitivity_to_accuracy(sensitivities)

  normal_noise_information <- data.frame()
  for (i in seq_along(sensitivities))
  {
    classifier <- get_normal_noise_classifier(sensitivity = sensitivities[i])
    row <- data.frame(accuracy    = accuracies[i],
                      sensitivity = sensitivities[i],
                      info        = get_information(classifier))
    normal_noise_information <- rbind(normal_noise_information, row)
  }

  normal_noise_information
}

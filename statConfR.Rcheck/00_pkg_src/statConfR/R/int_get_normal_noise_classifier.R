get_normal_noise_classifier <- function(sensitivity = NULL,
                                       accuracy    = NULL){
  if (is.null(sensitivity))
    sensitivity <- transform_normal_accuracy_to_sensitivity(accuracy)

  if (is.null(accuracy))
    accuracy <- transform_normal_sensitivity_to_accuracy(sensitivity)

  prior <- c(.5, .5)
  if (sensitivity == 0)         return(as.matrix(prior))
  if (is.infinite(sensitivity)) return(diag(prior))

  x  <- seq(-7 - sensitivity, 7 + sensitivity, length.out = 5000)
  delta <- x[2] - x[1]
  f1 <- dnorm(x, mean = -sensitivity/2) * delta * 1/2
  f2 <- dnorm(x, mean = +sensitivity/2) * delta * 1/2

  classifier <- rbind(f1, f2)
  classifier
}

transform_normal_accuracy_to_sensitivity <- function(accuracy){
  sensitivity <- 2 * qnorm(accuracy)
  sensitivity
}

transform_normal_sensitivity_to_accuracy <- function(sensitivity){
  accuracy <- pnorm(sensitivity / 2)
  accuracy
}



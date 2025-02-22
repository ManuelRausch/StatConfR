#' @importFrom utils tail
get_analytic_information_bounds <- function(prior, accuracies){
  upper <- sapply(accuracies, \(a) { get_upper_info_for_one(prior, a) })
  lower <- sapply(accuracies, \(a) { get_lower_info_for_one(prior, a) })

  information_bounds <- data.frame(accuracies = accuracies,
                                   highest    = upper     ,
                                   lowest     = lower     )
  information_bounds
}

get_upper_info_for_one <- function(prior, accuracy) {
  m1 <- floor(1/accuracy)
  m2 <- floor(1/accuracy)+1

  HY  <- get_entropy(prior)
  HYC <- ( (1/m1 - accuracy) * log(m2, 2) +
           (accuracy - 1/m2) * log(m1, 2)   ) / (1/m1 - 1/m2)

  info <- HY - HYC
  info
}

get_lower_info_for_one <- function(prior, accuracy) {
  a <- accuracy
  p <- sort(prior, decreasing = TRUE)
  L <- length(prior)

  # Edge case: If accuracy is 1, information will be equal to the overall
  # entropy because all uncertainty is reduced.
  if (a == 1)
  {
    info <- sum(p * log(1/p, 2))
    return(info)
  }

  # Edge case: If accuracy is as low as it can get (when simply predicting the
  # most probable label), information is zero.
  if (a == max(p))
  {
    info <- 0
    return(info)
  }

  # Determine which labels Y occur frequently enough (1:m3)
  s  <- p >= (cumsum(p) - a)/(1:L - 1)
  m3 <- tail(which(s), 1)
  if (length(m3) == 0) return(NaN)
  q <- sum(p[1:m3])
  pl <- p[1:m3]

  # Evaluate information for these labels because the others are ignored
  HY  <- sum(pl * log(1/pl, 2))
  HYC <- a * log(1/a, 2) + (q - a) * log((m3 - 1)/(q - a), 2)

  info <- HY - HYC
  names(info) <- NULL # Remove name inherited from m3
  info
}

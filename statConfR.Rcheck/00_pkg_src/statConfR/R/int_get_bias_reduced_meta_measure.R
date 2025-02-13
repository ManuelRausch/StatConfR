# Because Meta-I measures are inherently biased, simulate data to estimate and
# then subtract this bias.

#' @importFrom stats rmultinom

get_bias_reduced_meta_I_measures <- function(estimated_table)
{
  meta_I_measures <- get_meta_I_measures(estimated_table) # Baseline estimate

  # Simulations based on the observed frequencies
  nsim                    <- 1000
  simulated_meta_measures <- data.frame()
  for (i in 1:nsim)
  {
    # Simulate one data row-wise
    simulated_table <- estimated_table*0
    for (j in 1:nrow(estimated_table))
    {
      n <- sum(estimated_table[j, ])
      simulated_table[j, ] <- rmultinom(1, n, estimated_table[j, ]/n)
    }

    # Skip a simulation if accuracy is 50% or 100%
    estimated_classifier <- simulated_table/sum(simulated_table)
    a <- get_accuracy(estimated_classifier)
    if (round(a - 1, 6) == 0) next;
    if (round(a - 0, 6) == 0) next;

    simulated_meta_measures <- rbind(simulated_meta_measures,
                                     get_meta_I_measures(simulated_table))

    # Loading bar
    cat(sprintf('|%s%s|\r',
      paste0(rep('=', round(i/nsim*20)), collapse = ''),
      paste0(rep(' ', 20-round(i/nsim*20)), collapse = '')))
    if (i == nsim) cat("\n")
  }

  # Reduce bias
  expected_meta_I_measures     <- colMeans(simulated_meta_measures, na.rm = TRUE)
  estimated_bias               <- meta_I_measures - expected_meta_I_measures
  bias_reduced_meta_I_measures <- meta_I_measures - estimated_bias

  bias_reduced_meta_I_measures
}

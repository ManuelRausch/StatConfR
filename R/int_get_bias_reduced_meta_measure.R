# Because Meta-I measures are inherently biased, simulate data to estimate and
# then subtract this bias.
get_bias_reduced_meta_measure <- function(estimated_classifier, ns, get_meta_measure)
{
  meta_I_measure <- get_meta_measure(estimated_classifier) # Baseline estimate

  # Simulations based on the observed frequencies
  nsim                    <- 1000
  simulated_meta_measures <- c()
  for (i in 1:nsim)
  {
    # Simulate one data row-wise
    counts <- estimated_classifier*0
    for (j in 1:nrow(estimated_classifier)) 
    {
      counts[j, ] <- rmultinom(1, ns[j], estimated_classifier[j, ])
    }
    simulated_classifier <- counts/sum(counts)

    # Skip a simulation if accuracy is 50% or 100%
    a <- get_accuracy(simulated_classifier)
    if (round(a - 1, 6) == 0) next;
    if (round(a - 0, 6) == 0) next;

    simulated_meta_measures[i] <- get_meta_measure(simulated_classifier)
  }

  # If simulations did not work, return the baseline estimate
  if (length(simulated_meta_measures) < 1)
  {
    simulated_meta_measures <- na.omit(simulated_meta_measures)
    simulated_meta_measures <- get_meta_measure(estimated_classifier)
  }
  
  # Reduce bias
  expected_meta_measure       <- mean(simulated_meta_measures, na.rm = TRUE)
  estimated_bias              <- meta_I_measure - expected_meta_measure
  bias_reduced_meta_I_measure <- meta_I_measure - estimated_bias

  bias_reduced_meta_I_measure
}

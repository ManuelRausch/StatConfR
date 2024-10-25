# Summarize the frequencies of stimulus x response (decision and confidence
# rating) combinations in a contingency matrix with one row per stimulus and
# one column per response. Columnwise responses are sorted from "certaintly
# the first stimulus" (left) to "certainly the second stimulus" (right).
estimate_contingency_table <- function(data)
{
  number_of_confidence_steps <- length(levels(data$rating))
  confidence_steps           <- seq_len(number_of_confidence_steps)
  response_steps             <- c(-rev(confidence_steps), confidence_steps)

  estimated_table <- matrix(0, nrow = 2, ncol = 2*length(confidence_steps))
  rownames(estimated_table) <- paste0("Stimulus = ", c("-1", "+1"))
  colnames(estimated_table) <- paste0("Response = ", response_steps)
  
  tab <- table(data$stimulus, data$rating, data$correct)

  for (stimulus in levels(data$stimulus)) { # First stimulus = 0, second stimulus = 1
    for (correct in 0:1) { # Correct = 1, incorrect = 0
      for (rating in levels(data$rating)) { # Number of confidence bin
        # Selector for the target observations
        s <- data$stimulus == stimulus &
             data$correct  == correct  &
             data$rating   == rating

        which_stimulus <- which(stimulus == levels(data$stimulus))

        true_label <- 2 * which_stimulus - 3 # Into +1/-1 format
        decision   <- 2 * correct - 1 # Into +1/-1 format
        response <- decision * true_label * which(rating == levels(data$rating))

        which_response <- which(response == response_steps)

        estimated_table[which_stimulus, which_response] <- sum(s)
      }
    }
  }

  estimated_table
}

estimate_classifier <- function(dat)
{
  counts_table         <- get_counts_table(dat)
  estimated_classifier <- counts_table / nrow(dat)
  estimated_classifier
}

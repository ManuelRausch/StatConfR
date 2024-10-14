get_counts_table <- function(dat)
{
  if (!is.factor(dat$r)) stop("Variable r must be a factor with ordered levels: First half indicating y=-1 and second half indicating y=+1.")
  counts_table <- table(dat$y, dat$r)
  counts_table
}

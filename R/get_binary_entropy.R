get_binary_entropy <- H2 <- function(p)
{
  binary_entropy <- sapply(p, function(p) { get_entropy(c(p, 1-p)) } )
  binary_entropy
}


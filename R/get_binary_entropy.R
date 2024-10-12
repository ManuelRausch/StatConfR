#' Get Binary Entropy
#'
#' Compute entropy of a binary random variable where one realization has
#' probability p. Which of the two realiziations of a binary variable doesn't
#' matter because of symmetry, H2(p) = H2(1-p).
#' 
#' @param p Probability of one realiziation
#'
#' @return Entropy of the binary random variable
get_binary_entropy <- H2 <- function(p) 
{ 
  binary_entropy <- sapply(p, function(p) { get_entropy(c(p, 1-p)) } )
  binary_entropy
}


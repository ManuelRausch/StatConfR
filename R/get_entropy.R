#' Get Entropy 
#'
#' Compute entropy of a discrete random variable with probability vector p.
#' 
#' @param p Probability vector, has to sum to 1
#' @param base Base for the logarithm. If base = 2, the entropy is measured in
#'  bit. If base = e, the entropy is measured in nat. Default is base = 2.
#'
#' @return Entropy H(p)
get_entropy <- H <- function(p, base = 2)
{
  stopifnot( round(sum(p) - 1, 6) == 0 )

  p <- p[p != 0] # By convention, 0*log(1/0) = 0 
                # because lim_{p->0} p*log(1/p) = 0

  entropy <- sum( p * log(1/p, base) )
  entropy
}

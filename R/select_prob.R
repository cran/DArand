
#' Probabilities to select a normalization set without DE-gene
#'
#' @param m, number of genes
#' @param k, normalization subset size
#' @param invariant, boolean, when TRUE, probability of selection is evaluated for invariant gene
#'
#' @return a vector of probabilities of having at least one differential expression used as an reference selected in the normalization subset for any number of differential expressions d in the gene collection.
#' @export
#'
#' @examples
#'
#' select_prob(500, 10, invariant=TRUE)
select_prob = function(m,k, invariant=TRUE) {
  ok = rep(0,m-k+1)
  ok[1] = 1 # this is d=0
  for (d in 1:(m-k)) ok[d+1] = ok[d]*(m-k-d)/(m-d)
  if (!invariant) ok = c(1,ok[1:(m-k)])
  return(ok)
}

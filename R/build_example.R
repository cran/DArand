
#' Simulation of gene expressions using independant negative binomials
#'
#' @param m, number of genes
#' @param m1, number of differentially expressed genes. In the expression matrix, m1 first columns contain differentially expressed genes.
#' @param n1, number of samples under the first condition. The first n1 rows in the expression matrix.
#' @param n2, number of samples under the second condition (default n2=n1)
#' @param fold, maximal fold change added to the first m1 genes. The fold decreases proportionally to `1/sqrt(1:m1)`.
#' @param mu0, mean relative expression
#' @param use.scales, if TRUE random scales are used, otherwise all scales are set to 1.
#' @param nb.size, number of successful trials in the negative binomial distribution. If nb.size is set to Inf (default), the Poisson model is used.
#'
#' @details The function generates a list, of which the first element `X` is a matrix of `n1+n2` and m dimension with simulated expressions under Poisson or Negative Binomial distribution. Lines `1:n1` correspond to the first condition (or sub-group) and lines `(n1+1):(n1+n2)` to the second one. Columns `1:m1` contain counts imitating differential expressions.
#' @details In the ideal situation there is no microscopical variability between samples and all scales (so-called scaling factors) would be the same.  To simulate examples corresponding to this perfect situation, use argument `use.scales=FALSE` which will set all scales to 1.  When `use.scales=TRUE`, scales are simulated under uniform distribution *Unif(0.25,4)*.
#' @details The fold is maximal for the first expression and decreases proportionally to `1/sqrt(1:m1)`. The smallest fold  `fold/sqrt(m1)` is set to the *m1*-th expression.
#'
#' @return A list with components
#' \describe{
#'   \item{`X`}{a two-dimensional array containing the expression table of n individuals in rows and m gene expressions in columns. }
#'   \item{`m1`}{number of differentially expressed genes (as in arguments).}
#'   \item{`n1`}{number of samples under the first condition (as in arguments).}
#'   \item{`n2`}{number of samples under the second condition (as in arguments).}
#'   \item{`fold`}{maximal fold change between the differentally expressed genes and invariant genes (as in arguments). }
#'   \item{`scales`}{vector of simulated scales.}
#'   \item{`mu0`}{mean relative expression (as in arguments).}
#' }
#'
#' @export
#'
#' @examples
#'
#' L = build_example(m=500,m1=25,n1=6,fold=20,mu0=100,use.scales=FALSE,nb.size=Inf)
build_example = function(m=500,m1,n1=6,n2=n1,fold=100,mu0=100,use.scales=FALSE,nb.size=Inf) {

  n = n1+n2
  if (length(nb.size)==1) nb.size = rep(nb.size,l=2)

  if (use.scales) {
    scales = stats::runif(n,min=0.25,max=4)
    scales = scales/mean(scales)
    print(scales)
  } else scales=rep(1,n)

  # mean expression levels
  M = matrix(mu0, nrow=n, ncol=m)
  if (m1>0) for (j in 1:m1) M[1:n1,j] = M[1:n1,j] + fold/sqrt(j) # m1 DE-genes with fold-change

  M = M*scales
  if (nb.size[1]==Inf & nb.size[2]==Inf) { # Poisson case
    X = matrix(sapply(M,function(lambda) stats::rpois(n=1,lambda)), ncol=m, byrow=FALSE)
  } else { # Negative binomial
    A = matrix(sapply(M[1:n1,,drop=F],function(lambda) stats::rnbinom(n=1,nb.size[1],nb.size[1]/(nb.size[1]+lambda))), ncol=m, byrow=FALSE)
    B = matrix(sapply(M[-(1:n1),,drop=F],function(lambda) stats::rnbinom(n=1,nb.size[2],nb.size[2]/(nb.size[2]+lambda))), ncol=m, byrow=FALSE)
    X = rbind(A,B)
  }
  list(X=X,m1=m1,n1=n1,n2=n2,fold=fold,scales=scales,mu0=mu0)
}

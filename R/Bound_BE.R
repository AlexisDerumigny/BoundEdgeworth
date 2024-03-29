
#' Compute a Berry-Esseen-type bound
#'
#' This function returns a valid value \mjseqn{\delta_n} for the bound
#' \mjtdeqn{\sup_{x \in R}
#' \left| \textrm{Prob}(S_n \leq x) - \Phi(x) \right|
#'   \leq \delta_n,
#' }{\sup_{x \in \mathbb{R}}
#' \left| \textrm{Prob}(S_n \leq x) - \Phi(x) \right|
#'   \leq \delta_n,
#' }{sup_{x \in \mathbb{R}} | Prob(S_n <= x) - \Phi(x) | <= \delta_n,
#' }
#'
#' where \mjseqn{X_1, \dots, X_n} be \mjseqn{n} independent centered variables,
#' and \mjseqn{S_n} be their normalized sum, in the sense that
#' \mjseqn{S_n := \sum_{i=1}^n X_i / \textrm{sd}(\sum_{i=1}^n X_i)}.
#' This bounds follows from the triangular inequality
#' and the bound on the difference between a cdf and its 1st-order Edgeworth Expansion.
#'
#' \loadmathjax
#'
#' Note that the variables \mjseqn{X_1, \dots, X_n} must be independent
#' but may have different distributions (if \code{setup$iid = FALSE}).
#'
#'
#' @inheritParams Bound_EE1
#'
#' @return A vector of the same size as \code{n} with values \mjseqn{\delta_n}
#' such that
#' \mjtdeqn{\sup_{x \in R}
#' \left| \textrm{Prob}(S_n \leq x) - \Phi(x) \right|
#'   \leq \delta_n.
#' }{\sup_{x \in \mathbb{R}}
#' \left| \textrm{Prob}(S_n \leq x) - \Phi(x) \right|
#'   \leq \delta_n.
#' }{sup_{x \in R} | Prob(S_n <= x) - \Phi(x) | <= \delta_n.
#' }
#'
#'
#' @references Derumigny A., Girard L., and Guyonvarch Y. (2023).
#' Explicit non-asymptotic bounds for the distance to the first-order Edgeworth expansion,
#' Sankhya A. \doi{10.1007/s13171-023-00320-y}
#' \href{https://arxiv.org/abs/2101.05780}{arxiv:2101.05780}.
#'
#' @seealso \code{\link{Bound_EE1}()} for a bound on the distance
#' to the first-order Edgeworth expansion.
#'
#' @examples
#' setup = list(continuity = FALSE, iid = FALSE, no_skewness = FALSE)
#' regularity = list(C0 = 1, p = 2, kappa = 0.99)
#'
#' computedBound_EE1 <- Bound_EE1(
#'   setup = setup, n = 150, K4 = 9,
#'   regularity = regularity, eps = 0.1 )
#'
#' computedBound_BE <- Bound_BE(
#'   setup = setup, n = 150, K4 = 9,
#'   regularity = regularity, eps = 0.1 )
#'
#' print(c(computedBound_EE1, computedBound_BE))
#'
#' @export
#'
Bound_BE <- function(
  setup = list(continuity = FALSE, iid = FALSE, no_skewness = FALSE),
  n,
  K4 = 9, K3 = NULL, lambda3 = NULL, K3tilde = NULL,
  regularity = list(C0 = 1, p = 2),
  eps = 0.1)
{

  ub_DeltanE <- Bound_EE1(
    setup = setup, n = n,
    K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde,
    regularity = regularity, eps = eps)

  if (setup$no_skewness){

    ub_DeltanB <- ub_DeltanE

  } else {
    # If skewness, bounds on lambda3n is required.
    # It can be supplied by the user or obtained from K3 (or K4).

    env <- environment(); Update_bounds_on_moments(env)

    ub_DeltanB <- ub_DeltanE +
      abs(lambda3) * stats::dnorm(0, mean = 0, sd = 1) / (6  * sqrt(n))
  }

  return(ub_DeltanB)
}

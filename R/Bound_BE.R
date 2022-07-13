
#' Produce a Berry-Esseen-type bound
#'
#' This bounds follows from the triangular inequality
#' and the bound on the difference
#' between a cdf and its 1st-order Edgeworth Expansion
#'
#' @inheritParams Bound_EE1
#'
#' @examples
#' setup = list(continuity = FALSE, iid = FALSE, no_skewness = FALSE)
#' regularity = list(C0 = 1, p = 2, kappa = 0.99)
#'
#' #computedBound_EE1 <- Bound_EE1(
#' #  setup = setup, n = 150, K4 = 9,
#' #  regularity = regularity, eps = 0.1 )
#'
#' #computedBound_BE <- Bound_BE(
#' #  setup = setup, n = 150, K4 = 9,
#' #  regularity = regularity, eps = 0.1 )
#'
#' #print(c(computedBound_EE1, computedBound_BE))
#'
#' @export
#'
Bound_BE <- function(
  setup = list(continuity = FALSE, iid = FALSE, no_skewness = FALSE),
  n, K4, K3 = NULL, lambda3 = NULL, K3tilde = NULL,
  regularity = list(C0 = 1, p = 2, kappa = 0.99),
  eps = 0.1)
{

  ub_DeltanE <- Bound_EE1(
    setup = setup,
    n = n, K4 = K4,
    K3 = K3, lambda3 = lambda3, K3tilde = K3tilde,
    regularity = regularity,
    eps = eps)

  if (setup$no_skewness){

    ub_DeltanB <- ub_DeltanE

  } else {

    #--------------------------------------------------------------------------#
    # If skewness, required bounds on lambda3n, can be given
    # or obtained from K3, obtained from K4
    #--------------------------------------------------------------------------#

    # K3, lambda3, K3tilde optional since an upper bound on K4 is enough to
    # upper bound them.

    # Bound (by K4) on K3 if its value (or bound on it) is no provided
    if (is.null(K3)){
      K3 <- K4^(3/4)
    }

    # Bound (by K3) on lambda3 if its value (or bound on it) is no provided
    if (is.null(lambda3)){
      lambda3 <- Value_cst_bound_lambda3_by_K3() * K3
    }

    ub_DeltanB <- ub_DeltanE +
      abs(lambda3) * (dnorm(0, mean = 0, sd = 1) / 6) / sqrt(n)

  }

  return(ub_DeltanB)
}


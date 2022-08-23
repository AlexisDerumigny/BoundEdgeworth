
#------------------------------------------------------------------------------#
# Main function                                                            #####
#------------------------------------------------------------------------------#

#' Function that compute integrated R_n^iid
#' denoted with a bar: \eqn{\bar{R}_n^{iid}(eps)} in the paper
#' and defined in section C.3.2.
#'
#' @examples
#' Rn_iid_integrated(eps = 0.1, noskewness = TRUE, n = 300, K4 = 9)
#' Rn_iid_integrated(eps = 0.1, noskewness = FALSE, n = 300, K4 = 9, lambda3 = 6)
#'
#' @noRd
#'
Rn_iid_integrated <- function(eps, noskewness, n, K4, lambda3 = NULL)
{

  bound_modulus_psi <- Value_cst_bound_modulus_psi()

  denom <- 2 * (1 - 3*eps)^2 * pi

  tildeA2n <- bound_modulus_psi / denom * K4 / (24 * n^2) *
    2^(3) * Standard_gamma(5)

  tildeA5n <- bound_modulus_psi / denom * K4^2 / (576 * n^3) *
    2^(4) * Standard_gamma(6)

  val_e_2n <- e_2n(
    eps = eps, noskewness = noskewness, n = n, K4 = K4, lambda3 = lambda3)

  val_P_2n <- P_2n(
    eps = eps, noskewness = noskewness, n = n, K4 = K4, lambda3 = lambda3)

  tildeA6n <- bound_modulus_psi / pi * val_e_2n / (8 * n^2) *
    (K4 / 12 +  1 / (4 * (1 - 3*eps)^2) + val_P_2n / (576 * (1 - 3*eps)^2)) *
    2^(4) * Standard_gamma(6)

  if (isTRUE(noskewness)) {

    tildeA1n <- tildeA3n <- tildeA4n <- tildeA7n <- 0

  } else {

    tildeA1n <- bound_modulus_psi / denom * abs(lambda3) / (6 * n^(3/2)) *
      2^(0.5) * Standard_gamma(2.5)

    tildeA3n <- bound_modulus_psi / denom * lambda3^2 / (36 * n^2) *
      2^(3) * Standard_gamma(5)

    tildeA4n <- bound_modulus_psi / denom * K4 * abs(lambda3) / (72 * n^(5/2)) *
      2^(3.5) * Standard_gamma(5.5)

    tildeA7n <- bound_modulus_psi / pi * val_e_2n * abs(lambda3) / (12 * n^(3/2)) *
      (K4 / 12 +  1 / (4 * (1 - 3*eps)^2) + val_P_2n / (576 * (1 - 3*eps)^2)) *
      2^(3.5) * Standard_gamma(5.5)

  }

  return(tildeA1n + tildeA2n + tildeA3n + tildeA4n + tildeA5n + tildeA6n + tildeA7n)
}

#------------------------------------------------------------------------------#
# Auxiliary functions                                                      #####
#------------------------------------------------------------------------------#

#' Function that computes P_{2,n}(eps)
#' defined in section C.2 of the paper.
#'
#' @examples
#' P_2n(eps = 0.1, noskewness = TRUE, n = 500, K4 = 9)
#' P_2n(eps = 0.1, noskewness = FALSE, n = 500, K4 = 9, lambda3 = 6)
#'
#' @noRd
#'
P_2n <- function(eps, noskewness, n, K4, lambda3 = NULL){

  val_P_2n_noskewness <-
    48 * eps * sqrt(K4/n) +
    4 * eps^2 * K4 / n

  if (isTRUE(noskewness)) {

    val_P_2n <- val_P_2n_noskewness

  } else {

    val_P_2n <- val_P_2n_noskewness +
      96 * sqrt(2*eps) * abs(lambda3) / (K4^(1/4) * n^(1/4)) +
      32 * eps * lambda3^2 / sqrt(K4 * n) +
      16 * sqrt(2) * K4^(1/4) * abs(lambda3) * eps^(3/2) / n^(3/4)
  }

  return (val_P_2n)
}

#' Function that computes e_{2,n}(eps)
#' defined in section C.2 of the paper.
#'
#' @examples
#' e_2n(eps = 0.1, noskewness = TRUE, n = 500, K4 = 9)
#' e_2n(eps = 0.1, noskewness = FALSE, n = 500, K4 = 9, lambda3 = 6)
#'
#' @noRd
#'
e_2n <- function(eps, noskewness, n, K4, lambda3 = NULL){

  val_P_2n <- P_2n(
    eps = eps, noskewness = noskewness, n = n, K4 = K4, lambda3 = lambda3)

  val_e_2n <- exp(eps^2 *
                    (1/6 + 1/(2 * (1 - 3*eps)^2) +
                       2 * val_P_2n / (576 * (1 - 3*eps)^2)))

  return(val_e_2n)
}

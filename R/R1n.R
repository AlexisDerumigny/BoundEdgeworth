

#' Function that computes Rn_inid_integrated(eps) (denoted with a bar)
#' Defined in C.3.1
#'
#' @examples
#' Rn_inid_integrated(eps = 0.1, noskewness = TRUE, n = 200, K4 = 9)
#' Rn_inid_integrated(eps = 0.1, noskewness = FALSE, n = 200, K4 = 9, lambda3 = 6)
#'
#' @noRd
#'
Rn_inid_integrated <- function(eps, noskewness, n, K4, lambda3 = NULL)
{

  bound_modulus_psi <- Value_cst_bound_modulus_psi()

  A1n <- bound_modulus_psi / (48 * (1 - 3 * eps)^2 * pi) * (K4 / n)^(3/2)
    2^4 * Standard_gamma(5)

  A2n <- bound_modulus_psi / (24^2 * 2 * (1 - 3 * eps)^2 * pi) * (K4 / n)^2 *
    2^(5) * Standard_gamma(6)

  A6n <- bound_modulus_psi * e_1n(eps = eps, noskewness = noskewness) / (2 * pi) *
    (K4 / n)^2 *
    (1/24 + P_1n(eps = eps, noskewness = noskewness) / (2 * (1 - 3 * eps)^2))^2 *
    2^5 * Standard_gamma(6)

  if (isTRUE(noskewness)) {

    A3n <- A4n <- A5n <- A7n <- 0

  } else {

    A3n <- bound_modulus_psi / (12 * (1 - 3 * eps)^2 * pi) * (K4 / n)^(5/4) *
      2^(3.5) * Standard_gamma(4.5)

    A4n <- bound_modulus_psi / (72 * (1 - 3 * eps)^2 * pi) * (K4 / n)^(3/2) *
      2^4 * Standard_gamma(5)

    A5n <- bound_modulus_psi / (144 * (1 - 3 * eps)^2 * pi) * (K4 / n)^(7/4) *
      2^(4.5) * Standard_gamma(5.5)

    A7n <- bound_modulus_psi * e_1n(eps = eps, noskewness = noskewness) / (6 * pi) *
      abs(lambda3) * K4 / n^(3/2) *
      (1/24 + P_1n(eps = eps, noskewness = noskewness) / (2 * (1 - 3 * eps)^2)) *
      2^5 * Standard_gamma(6)

  }

  result = A1n + A2n + A3n + A4n + A5n + A6n + A7n

  return (result)
}

#' Function that computes P_{1,n}(eps)
#' Defined in C.3.1
#'
#' @examples
#' P_1n(eps = 0.1, noskewness = TRUE)
#' P_1n(eps = 0.1, noskewness = FALSE)
#'
#' @noRd
#'
P_1n <- function(eps, noskewness){

  val_P_1n_noskewness <-
    (144 + 48*eps + 4*eps^2) / 576

  if (isTRUE(noskewness)){

    val_P_1n <- val_P_1n_noskewness

  } else {

    additional_part_skewness <-
      (96*sqrt(2*eps) + 32*eps + 16*sqrt(2)*eps^(3/2)) / 576

    val_P_1n <- val_P_1n_noskewness + additional_part_skewness
  }

  return(val_P_1n)
}

#' Function that computes e_{1,n}(eps)
#' Defined in C.3.1
#'
#' @examples
#' e_1n(eps = 0.1, noskewness = TRUE)
#' e_1n(eps = 0.1, noskewness = FALSE)
#'
#' @noRd
#'
e_1n <- function(eps, noskewness){

  val_P_1n <- P_1n(eps = eps, noskewness = noskewness)

  val_e_1n <- exp(eps^2 * (1/6 + 2*val_P_1n / (1 - 3*eps)^2))

  return(val_e_1n)
}

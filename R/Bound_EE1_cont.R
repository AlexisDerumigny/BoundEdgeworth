


# Main functions -------------------------------------------------

Bound_EE1_cont_inid_skew <- function(n, eps, K4, K3 = NULL, lambda3 = NULL, K3tilde = NULL)
{

}

Bound_EE1_cont_inid_noskew <- function(n, eps, K4, K3 = NULL, lambda3 = NULL, K3tilde = NULL)
{

}

Bound_EE1_cont_iid_skew <- function(n, eps, K4, K3 = NULL, lambda3 = NULL, K3tilde = NULL)
{

}

Bound_EE1_cont_iid_noskew <- function(n, eps, K4, K3 = NULL, lambda3 = NULL, K3tilde = NULL)
{

}


# Remainder terms  ------------------------------------------------------



#------------------------------------------------------------------------------#
##### Continuity, inid - Theorem A.4 (Thm 3.1 under Assumption 2.1) #####
#------------------------------------------------------------------------------#

#' Compute the curly brace related to Delta for
#' r_{1,n}^{inid, skew} and r_{1,n}^{inid, noskew}
#'
#' @param p : parameter not introduced in the paper but to wrap-up
#' the two expressions for r_{1,n}^{inid, skew} and r_{1,n}^{inid, noskew}
#' p = 3 for r_{2,n}^{inid, skew}
#' p = 4 for r_{1,n}^{inid, noskew}
#'
#' @noRd
#'
Delta_curly_brace_part_r2n <- function(eps, p, n, K4, K3tilde){

  Delta <- Delta_of_K4_and_n(K4 = K4, n = n)

  if (Delta == 0) {

    upper_end <- 16 * pi^3 * n^2 / K3tilde^4

    lower_end_if_integral_not_nul <- sqrt(2*eps) * (n/K4)^(1/4)

    if (upper_end <= lower_end_if_integral_not_nul) {
      value <- 0
    } else {
      value <- (upper_end^p - lower_end_if_integral_not_nul^p) / p
    }

  } else {

    value <- 0.5 * abs(Delta)^(-1 - p/2) * abs(
      Lower_Incomplete_gamma(p / 2, 2^8 * pi^6 * Delta * n^4 / K3tilde^8) -
        Lower_Incomplete_gamma(p / 2,
                               Delta *
                                 min(2 * eps * sqrt(n/K4),
                                     2^8 * pi^6 * n^4 / K3tilde^8)))
  }

  return(value)
}

r2n_inid_skew <- function(eps, n, lambda3, K3tilde, K4)
{
  value_Rn_inid_integrated <- Rn_inid_integrated(
    eps = eps, n = n, K4 = K4, lambda3 = lambda3, noskewness = FALSE)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r2n(
    eps = eps, p = 3, n = n, K4 = K4, K3tilde = K3tilde)

  min_lower_end_Gamma <- min(sqrt(2*eps) * (n/K4)^(1/4), 16 *pi^3 * n^2 / K3tilde^4)

  bound_modulus_psi <- Value_cst_bound_modulus_psi()

  t1star <- Value_t1star()

  shortcut <- 2 * n * (1 - 4 * pi * Value_chi1() * t1star) / K3tilde^2

  value_r2n_inid_skew <-
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
    0.3334 * K3tilde^4 * abs(lambda3) / (16 * pi^4 * n^(5/2)) +
    14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
    4.3394 * K3tilde^12 * abs(lambda3) / ((2*pi)^12 * n^(13/2)) +
    abs(lambda3) *
    (Upper_Incomplete_gamma(3/2, min_lower_end_Gamma) -
       Upper_Incomplete_gamma(3/2, 16 *pi^3 * n^2 / K3tilde^4)) / sqrt(n) +
    value_Rn_inid_integrated +
    bound_modulus_psi * K3 * value_Delta_curly_brace_part / (6 * pi * sqrt(n)) +
    bound_modulus_psi / (2*pi) *
    (Standard_gamma(shortcut) - Standard_gamma(shortcut * (pi * t1star)^2)) +
    bound_modulus_psi / (4*pi) *
    abs(Upper_Incomplete_gamma(0, (2*pi)^7 * n^4 / K3tilde^8) -
          Upper_Incomplete_gamma(0, 2 * n / K3tilde^2))
}

r2n_inid_noskew <- function(eps, n, lambda3, K3tilde, K4)
{
  value_Rn_inid_integrated <- Rn_inid_integrated(
    eps = eps, noskewness = TRUE, n = n, K4 = K4)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r2n(
    eps = eps, p = 4, n = n, K4 = K4, K3tilde = K3tilde)

  min_lower_end_Gamma <- min(sqrt(2*eps) * (n/K4)^(1/4), 16 *pi^3 * n^2 / K3tilde^4)

  bound_modulus_psi <- Value_cst_bound_modulus_psi()

  t1star <- Value_t1star()

  shortcut <- 2 * n * (1 - 4 * pi * Value_chi1() * t1star) / K3tilde^2

  value_r2n_inid_noskew <-
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
    14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
    value_Rn_inid_integrated +
    bound_modulus_psi * K4 * value_Delta_curly_brace_part / (6 * pi * n) +
    bound_modulus_psi / (2*pi) *
    (Standard_gamma(shortcut) - Standard_gamma(shortcut * (pi * t1star)^2)) +
    bound_modulus_psi / (4*pi) *
    abs(Upper_Incomplete_gamma(0, (2*pi)^7 * n^4 / K3tilde^8) -
          Upper_Incomplete_gamma(0, 2 * n / K3tilde^2))
}






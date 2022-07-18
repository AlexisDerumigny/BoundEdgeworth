
# These are the functions to compute the bounds
# for Theorem A.4 and A.5

# Structure of this file:
# - main functions
# - functions for Theorem A.2
# - functions for Theorem A.3


#------------------------------------------------------------------------------#
#      Main functions                                                      #####
#------------------------------------------------------------------------------#

#' @examples
#' Bound_EE1_cont_inid_skew(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, lambda3 = 1, K3 = 4)
#'
#' @noRd
Bound_EE1_cont_inid_skew <- function(n, eps, K4, K3, lambda3, K3tilde)
{
  dominant_term <-
    Bound_EE1_cont_common_part_noskewness(eps = eps, n = n,
                                          K4 = K4, K3tilde = K3tilde)
  additional_term_skewness <-
    Bound_EE1_cont_additional_skewness(eps = eps, n = n,
                                       K3tilde = K3tilde, lambda3 = lambda3)

  remainder_term <-
    r2n_inid_skew(eps = eps, n = n, lambda3 = lambda3,
                  K3tilde = K3tilde, K4 = K4, K3 = K3)

  return(dominant_term + additional_term_skewness + remainder_term)
}

#' @examples
#' Bound_EE1_cont_inid_noskew(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, K3 = 4)
#'
#' @noRd
Bound_EE1_cont_inid_noskew <- function(n, eps, K4, K3, K3tilde)
{
  dominant_term <-
    Bound_EE1_cont_common_part_noskewness(eps = eps, n = n,
                                          K4 = K4, K3tilde = K3tilde)

  remainder_term <-
    r2n_inid_noskew(eps = eps, n = n,
                    K3tilde = K3tilde, K4 = K4)

  return(dominant_term + remainder_term)
}

#' @examples
#' Bound_EE1_cont_iid_skew(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, lambda3 = 1, K3 = 4)
#'
#' @noRd
Bound_EE1_cont_iid_skew <- function(n, eps, K4, K3, lambda3, K3tilde)
{
  dominant_term <-
    Bound_EE1_cont_common_part_noskewness(eps = eps, n = n,
                                          K4 = K4, K3tilde = K3tilde)
  additional_term_skewness <-
    Bound_EE1_cont_additional_skewness(eps = eps, n = n,
                                       K3tilde = K3tilde, lambda3 = lambda3)

  remainder_term <-
    r2n_iid_skew(eps = eps, n = n, lambda3 = lambda3,
                 K3tilde = K3tilde, K4 = K4, K3 = K3)

  return(dominant_term + additional_term_skewness + remainder_term)
}

#' @examples
#' Bound_EE1_cont_iid_noskew(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, K3 = 4)
#'
#' @noRd
Bound_EE1_cont_iid_noskew <- function(n, eps, K4, K3, K3tilde)
{
  dominant_term <-
    Bound_EE1_cont_common_part_noskewness(eps = eps, n = n,
                                          K4 = K4, K3tilde = K3tilde)

  remainder_term <-
    r2n_iid_noskew(eps = eps, n = n,
                   K3tilde = K3tilde, K4 = K4, K3 = K3)

  return(dominant_term + remainder_term)
}


#------------------------------------------------------------------------------#
#      Main terms                                                          #####
#------------------------------------------------------------------------------#

#' Main term of Theorem A.4 and A.5 (noskewness case)
#'
#' @examples
#' Bound_EE1_cont_common_part_noskewness(n = 300, eps = 0.1, K4 = 9, K3tilde = 6)
#'
#' @noRd
Bound_EE1_cont_common_part_noskewness <- function(n, eps, K4, K3tilde)
{
  return(
    0.327 * K4 * (1/12 + 1 / (4 * (1 - 3*eps)^2)) / n
  )
}

#' Additional main term for Theorem A.2 and A.3 in case of skewness
#'
#' @examples
#' Bound_EE1_cont_additional_skewness(n = 300, eps = 0.1, K3tilde = 6, lambda3 = 0.1)
#'
#' @noRd
Bound_EE1_cont_additional_skewness <- function(n, eps, K3tilde, lambda3)
{
  return(
    0.037 * e_1n(eps = eps, noskewness = FALSE) * lambda3^2 / n
  )
}

#------------------------------------------------------------------------------#
#      Remainder terms                                                     #####
#------------------------------------------------------------------------------#


#' @examples
#' r2n_inid_skew(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, lambda3 = 1, K3 = 4)
#'
#' @noRd
r2n_inid_skew <- function(n, eps, K4, K3, lambda3, K3tilde)
{
  value_Rn_inid_integrated <- Rn_inid_integrated(
    eps = eps, n = n, K4 = K4, lambda3 = lambda3, noskewness = FALSE)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r2n(
    eps = eps, p = 3, n = n, K4 = K4, K3tilde = K3tilde)

  upper_end_Gamma <- 16 *pi^3 * n^2 / K3tilde^4
  lower_end_Gamma <- min( sqrt(2*eps) * (n/K4)^(1/4) , upper_end_Gamma )

  result <-
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
    0.3334 * K3tilde^4 * abs(lambda3) / (16 * pi^4 * n^(5/2)) +
    14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
    4.3394 * K3tilde^12 * abs(lambda3) / ((2*pi)^12 * n^(13/2)) +
    abs(lambda3) * abs(Upper_Incomplete_gamma(3/2, upper_end_Gamma) -
                         Upper_Incomplete_gamma(3/2, lower_end_Gamma) ) / sqrt(n) +
    value_Rn_inid_integrated +
    Value_cst_bound_modulus_psi() * K3 * value_Delta_curly_brace_part / (6 * pi * sqrt(n)) +
    integral_terms_r2n(n = n, K3tilde = K3tilde)

  return (result)
}

r2n_inid_noskew <- function(eps, n, K3tilde, K4)
{
  value_Rn_inid_integrated <- Rn_inid_integrated(
    eps = eps, noskewness = TRUE, n = n, K4 = K4)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r2n(
    eps = eps, p = 4, n = n, K4 = K4, K3tilde = K3tilde)

  result <-
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
    14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
    value_Rn_inid_integrated +
    Value_cst_bound_modulus_psi() * K4 * value_Delta_curly_brace_part / (6 * pi * n) +
    integral_terms_r2n(n = n, K3tilde = K3tilde)

  return (result)
}


r2n_iid_skew <- function(eps, n, lambda3, K3tilde, K4, K3)
{
  value_Rn_iid_integrated <- Rn_iid_integrated(
    eps = eps, n = n, K4 = K4, lambda3 = lambda3, noskewness = FALSE)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r2n(
    eps = eps, p = 3, n = n, K4 = K4, K3tilde = K3tilde)

  upper_end_Gamma <- 16 *pi^3 * n^2 / K3tilde^4
  lower_end_Gamma <- min( sqrt(2*eps) * (n/K4)^(1/4) , upper_end_Gamma )

  result <-
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
    0.3334 * K3tilde^4 * abs(lambda3) / (16 * pi^4 * n^(5/2)) +
    14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
    4.3394 * K3tilde^12 * abs(lambda3) / ((2*pi)^12 * n^(13/2)) +
    abs(lambda3) * abs(Upper_Incomplete_gamma(3/2, upper_end_Gamma) -
                         Upper_Incomplete_gamma(3/2, lower_end_Gamma) ) / sqrt(n) +
    value_Rn_iid_integrated +
    integral_terms_r2n(n = n,
                       K3tilde = K3tilde) +
    1.306 * ( e_2n(eps = eps, noskewness = FALSE , n = n, K4 = K4, lambda3 = lambda3)
              - e3(eps = eps) ) * lambda3^2 / (36 * n)

  return (result)
}


r2n_iid_noskew <- function(eps, n, K3tilde, K4, K3)
{
  value_Rn_iid_integrated <- Rn_iid_integrated(
    eps = eps, n = n, K4 = K4, noskewness = TRUE)

  upper_end_Gamma <- 2^5 * pi^6 * n^4 / K3^8
  lower_end_Gamma <- min( eps * sqrt(n / (16 * K4)) , upper_end_Gamma)

  result <-
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
    14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
    value_Rn_iid_integrated +
    16 * 1.0253 * K3 * abs( Upper_Incomplete_gamma(2, upper_end_Gamma) -
                              Upper_Incomplete_gamma(2, lower_end_Gamma)
    ) /(3 * pi * n) +
    integral_terms_r2n(n = n, K3tilde = K3tilde)

  return (result)
}


#------------------------------------------------------------------------------#
#   Helper functions for the remainder terms                               #####
#------------------------------------------------------------------------------#

#' Compute the curly brace related to Delta for
#' r_{2,n}^{inid, skew} and r_{2,n}^{inid, noskew}
#'
#' @param p : parameter not introduced in the paper but to wrap-up
#' the two expressions for r_{2,n}^{inid, skew} and r_{2,n}^{inid, noskew}
#' p = 3 for r_{2,n}^{inid, skew}
#' p = 4 for r_{2,n}^{inid, noskew}
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


#' @examples
#' integral_terms_r2n(n = 300, K3tilde = 6)
#'
#' @noRd
integral_terms_r2n <- function(n, K3tilde)
{
  bound_modulus_psi <- Value_cst_bound_modulus_psi()

  t1star <- Value_t1star()

  shortcut <- 2 * n * (1 - 4 * pi * Value_chi1() * t1star) / K3tilde^2

  first_integral <- bound_modulus_psi / (2*pi) *
    abs(Upper_Incomplete_gamma(0, shortcut) -
          Upper_Incomplete_gamma(0, shortcut * (pi * t1star)^2))

  second_integral <- bound_modulus_psi / (2*pi) *
    abs(Upper_Incomplete_gamma(0, (2*pi)^7 * n^4 / K3tilde^8) -
          Upper_Incomplete_gamma(0, 2 * n / K3tilde^2))

  return (first_integral + second_integral)
}


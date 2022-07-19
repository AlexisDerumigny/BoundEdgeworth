
# These are the functions to compute the bounds
# for Theorem A.2 and A.3

# Structure of this file:
# - main terms of the expansion
# - functions for Theorem A.2
# - functions for Theorem A.3


#------------------------------------------------------------------------------#
#      Main terms                                                          #####
#------------------------------------------------------------------------------#

# Main term of Theorem A.2 and A.3, excluding skewness
Bound_EE1_nocont_common_part_noskewness <- function(eps, n, K4, K3tilde)
{
  return(
    0.1995 * K3tilde / sqrt(n) +
      (0.031 * K3tilde^2 + 0.327 * K4 * (1/12 + 1 / (4 * (1 - 3*eps)^2))) / n
  )
}

# Additional main term for Theorem A.2 and A.3 in case of skewness
Bound_EE1_nocont_additional_skewness <- function(eps, n, K3tilde, lambda3)
{
  return((0.054 * abs(lambda3) * K3tilde +
            0.037 * e_1n(eps = eps, noskewness = FALSE) * lambda3^2) / n
  )
}



#------------------------------------------------------------------------------#
#      Nocontinuity, inid - Theorem A.2 (Thm 2.1 under Assumption 2.1)     #####
#------------------------------------------------------------------------------#

## Main functions -------------------------------------------------

#' @examples
#' Bound_EE1_nocont_inid_skew(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, lambda3 = 1, K3 = 4)
#'
#' @noRd
Bound_EE1_nocont_inid_skew <- function(n, eps, K4, K3, K3tilde, lambda3)
{
  dominant_term <-
    Bound_EE1_nocont_common_part_noskewness(eps = eps, n = n,
                                            K4 = K4, K3tilde = K3tilde)
  additional_term_skewness <-
    Bound_EE1_nocont_additional_skewness(eps = eps, n = n,
                                         K3tilde = K3tilde, lambda3 = lambda3)

  remainder_term <-
    r1n_inid_skew(eps = eps, n = n, lambda3 = lambda3,
                  K3tilde = K3tilde, K4 = K4, K3 = K3)

  return(dominant_term + additional_term_skewness + remainder_term)
}

#' @examples
#' Bound_EE1_nocont_inid_noskew(n = 300, eps = 0.1, K4 = 9, K3tilde = 6)
#'
#' @noRd
Bound_EE1_nocont_inid_noskew <- function(n, eps, K4, K3tilde)
{
  dominant_term <-
    Bound_EE1_nocont_common_part_noskewness(eps = eps, n = n, K4 = K4, K3tilde = K3tilde)

  remainder_term <-
    r1n_inid_noskew(eps = eps, n = n, K4 = K4, K3tilde = K3tilde)

  return(dominant_term + remainder_term)
}


## Remainder terms -----------------------------------------------

r1n_inid_skew <- function(eps, n, lambda3, K3tilde, K4, K3)
{
  value_Rn_inid_integrated <- Rn_inid_integrated(
    eps = eps, n = n, K4 = K4, lambda3 = lambda3, noskewness = FALSE)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r1n(
    eps = eps, p = 3, n = n, K4 = K4, K3tilde = K3tilde)

  min_lower_end_Gamma <- min(sqrt(2*eps) * (n/K4)^(1/4), 2*sqrt(n) / K3tilde)

  value_r1n_inid_skew <-
    (14.1961 + 67.0415) * K3tilde^4 / (16 * pi^4 * n^2) +
    4.3394 * abs(lambda3) * K3tilde^3 / (8 * pi^3 * n^2) +
    value_Rn_inid_integrated +
    abs(lambda3) *
    (Upper_Incomplete_gamma(3/2, min_lower_end_Gamma) -
       Upper_Incomplete_gamma(3/2, 2*sqrt(n) / K3tilde)) / sqrt(n) +
    Value_cst_bound_modulus_psi() * K3 * value_Delta_curly_brace_part / (6 * pi * sqrt(n))
}

r1n_inid_noskew <- function(eps, n, K4, K3tilde)
{
  value_Rn_inid_integrated <- Rn_inid_integrated(
    eps = eps, noskewness = TRUE, n = n, K4 = K4)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r1n(
    eps = eps, p = 4, n = n, K4 = K4, K3tilde = K3tilde)

  min_lower_end_Gamma <- min(sqrt(2*eps) * (n/K4)^(1/4), 2*sqrt(n) / K3tilde)

  value_r1n_inid_skew <-
    (14.1961 + 67.0415) * K3tilde^4 / (16 * pi^4 * n^2) +
    value_Rn_inid_integrated +
    Value_cst_bound_modulus_psi() * K4 * value_Delta_curly_brace_part / (6 * pi * n)
}



#------------------------------------------------------------------------------#
#       Nocontinuity, iid - Theorem A.3 (Thm 2.1 under Assumption 2.2)     #####
#------------------------------------------------------------------------------#


## Main functions -----------------------------------------------------

Bound_EE1_nocont_iid_skew <- function(n, eps, K4, K3, lambda3, K3tilde)
{
  dominant_term <-
    Bound_EE1_nocont_common_part_noskewness(eps = eps, n = n,
                                            K4 = K4, K3tilde = K3tilde)
  additional_term_skewness <-
    Bound_EE1_nocont_additional_skewness(eps = eps, n = n,
                                         K3tilde = K3tilde, lambda3 = lambda3)

  remainder_term <-
    r1n_iid_skew(eps = eps, n = n, lambda3 = lambda3,
                 K3tilde = K3tilde, K4 = K4, K3 = K3)

  return (dominant_term + additional_term_skewness + remainder_term)
}

Bound_EE1_nocont_iid_noskew <- function(n, eps, K4, K3tilde)
{
  dominant_term <-
    Bound_EE1_nocont_common_part_noskewness(eps = eps, n = n, K4 = K4, K3tilde = K3tilde)

  remainder_term <- r1n_iid_noskew(eps = eps, n = n, K4 = K4, K3tilde = K3tilde)

  return(dominant_term + remainder_term)
}


## Remainder terms -------------------------------------------------------


#' Function to compute r_{1,n}^{iid, skew}(eps) defined in Section A.3
#'
#' @example
#' r1n_iid_skew(eps = 0.1, n = 200, K4 = 9, K3 = 5, K3tilde = 7, lambda3 = 6)
#'
#' @noRd
#'
r1n_iid_skew <- function(eps, n, K4, K3, K3tilde, lambda3)
{
  value_Rn_iid_integrated <- Rn_iid_integrated(
    eps = eps, noskewness = FALSE, n = n, K4 = K4, lambda3 = lambda3)

  value_e2n <- e_2n(
    eps = eps, noskewness = FALSE, n = n, K4 = K4, lambda3 = lambda3)

  min_lower_end_Gamma <- min(sqrt(2*eps) * (n/K4)^(1/4), 2*sqrt(n) / K3tilde)

  value_r1n_iid_skew <-
    (14.1961 + 67.0415) * K3tilde^4 / (16 * pi^4 * n^2) +
    4.3394 * abs(lambda3) * K3tilde^3 / (8 * pi^3 * n^2) +
    value_Rn_iid_integrated +
    1.306 * (value_e2n - e3(eps)) * lambda3^2 / (36 * n) +
    abs(lambda3) *
    (Upper_Incomplete_gamma(3/2, min_lower_end_Gamma) -
       Upper_Incomplete_gamma(3/2, 2*sqrt(n) / K3tilde)) / sqrt(n) +
    Value_cst_bound_modulus_psi() * 2^(5/2) * K3 / (3 * pi * sqrt(n)) *
    (Upper_Incomplete_gamma(3/2, 4*n / (8 * K3tilde^2)) -
       Upper_Incomplete_gamma(3/2, min_lower_end_Gamma^2 / 8))

  return(value_r1n_iid_skew)
}

#' Function to compute r_{1,n}^{iid, noskew}(eps) defined in Section A.3
#'
#' @examples
#' r1n_iid_noskew(eps = 0.1, n = 300, K4 = 9, K3tilde = 6)
#' TODO: negatif pour cet exemple; probleme quelque part?
#'
#' @noRd
#'
r1n_iid_noskew <- function(eps, n, K4, K3tilde)
{
  value_Rn_iid_integrated <- Rn_iid_integrated(
    eps = eps, noskewness = TRUE, n = n, K4 = K4)

  min_lower_end_Gamma <- min(sqrt(2*eps) * (n/K4)^(1/4), 2*sqrt(n) / K3tilde)

  value_r1n_iid_noskew <-
    (14.1961 + 67.0415) * K3tilde^4 / (16 * pi^4 * n^2) +
    value_Rn_iid_integrated +
    16 * Value_cst_bound_modulus_psi() * K4 / (3 * pi * n) *
    (Upper_Incomplete_gamma(2, 4*n / (8 * K3tilde^2)) -
       Upper_Incomplete_gamma(2, min_lower_end_Gamma^2 / 8))

  return(value_r1n_iid_noskew)
}



#------------------------------------------------------------------------------#
#####  #####
#------------------------------------------------------------------------------#



#' Compute the quantity \Delta that intervenes in the expression
#' of r_{1,n}^{inid, skew} and r_{1,n}^{inid, noskew}
#'
#' @noRd
Delta_of_K4_and_n <- function(K4, n){

  Delta <- (1 - 4*Value_chi1() - sqrt(K4/n)) / 2

  return(Delta)
}

#' Compute the curly brace related to Delta for
#' r_{1,n}^{inid, skew} and r_{1,n}^{inid, noskew}
#'
#' @param p : parameter not introduced in the paper but to wrap-up
#' the two expressions for r_{1,n}^{inid, skew} and r_{1,n}^{inid, noskew}
#' p = 3 for r_{1,n}^{inid, skew}
#' p = 4 for r_{1,n}^{inid, noskew}
#'
#' @noRd
#'
Delta_curly_brace_part_r1n <- function(eps, p, n, K4, K3tilde){

  Delta <- Delta_of_K4_and_n(K4 = K4, n = n)

  if (Delta == 0) {

    upper_end <- 2 * sqrt(n) / K3tilde

    lower_end_if_integral_not_nul <- sqrt(2*eps) * (n/K4)^(1/4)

    if (upper_end <= lower_end_if_integral_not_nul) {
      value <- 0
    } else {
      value <- (upper_end^p - lower_end_if_integral_not_nul^p) / p
    }

  } else {
    upper_end = 2 * Delta * min(eps * sqrt(n/K4), 2 * n / K3tilde^2)
    lower_end = 4 * Delta * n / K3tilde^2

    value <- 0.5 * abs(Delta)^(- p/2) *
      abs( Lower_Incomplete_gamma(p / 2, lower_end) -
             Lower_Incomplete_gamma(p / 2, upper_end) )
  }

  return(value)
}



#' Function to computes e_3(eps)
#' defined in Section A.3
#'
#' @example
#' e3(eps = 0.1)
#'
#' @noRd
#'
e3 <- function(eps)
{
  return(exp(eps^2 / 6 + eps^2 / (2 * (1 - 3*eps))^2))
}



# # Old code ?
#
# r1n_iid_skew <- function(eps, n, lambda3, K3tilde, K4)
# {
#   Rn_inid_integrated_of_eps <- Rn_inid_integrated(
#     eps = eps, n = n, K4 = K4, lambda3 = lambda3, noskewness = FALSE)
#
#
#
#
#   stop("replace the p by its value")
#   stop("fix the \\Gamma and \\gamma")
#
#   min = min( sqrt(2 * eps) * (n/K4)^(1/4) , 2 * sqrt(n) / K3tilde )
#
#   rninidskew = (14.1961 + 67.0415) * K3tilde^4 / (16 * pi^4 * n^2) +
#     4.3394 * abs(lambda3) * K3tilde^3 / (8 * pi^3 * n^2) +
#     R1int +
#     abs(lambda3) *
#     ( Gamma( 3/2 , min )
#       - Gamma( 3/2 , 2 * sqrt(n) / K3tilde) )  / sqrt(n) +
#
#     1.0253 * K3 / (6 * pi * sqrt(n)) *
#     ( 0.5 * abs(Delta)^(-1-3/2) *
#         abs( gamma(3/2, 4 * Delta * n / K3tilde^2)
#              - gamma(3/2, Delta * min^2 ) ) )
#
#   result = 1.2533 * K3tilde / (2 * pi * sqrt(n)) +
#     0.3334 * abs(lambda3) * K3tilde / (2 * pi * n) +
#     1.2187 * K3tilde^2 / (4 * pi * n) +
#     (0.327 * K4 / n) * (1/12 + 1 / (4 * (1 - 3 * eps)^2 ) ) +
#     (1.306 * e_1(eps = eps, noskewness = FALSE) * lambda3^2) / (36 * n) +
#     rninidskew
# }



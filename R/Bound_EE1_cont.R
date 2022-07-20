
# These are the functions to compute the bounds for Theorem A.4 and A.5
# WITHOUT the term involving the integral of f_{S_n}
# 1.0253 / pi * int_{a_n}^{b_n} |f_{S_n}(t)| / t dt
# Hence the suffix: _wo_int_fSn for these four main functions
# This term is dealt with given the regularity conditions in the final
# user function Bound_EE1().

#------------------------------------------------------------------------------#
#      Main functions                                                      #####
#------------------------------------------------------------------------------#

#' Compute the bound of Theorem A.4
#' (Theorem 3.1 under Assumption 2.1: inid case)
#' in the skewness case (Equation (23))
#' WITHOUT the term involving the integral of f_{S_n}
#'
#' @examples
#' Bound_EE1_cont_inid_skew_wo_int_fSn(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, lambda3 = 1, K3 = 4)
#'
#' @noRd
#'
Bound_EE1_cont_inid_skew_wo_int_fSn <- function(n, eps, K4, K3, lambda3, K3tilde)
{
  main_term <-
    Bound_EE1_cont_common_part_noskewness(eps = eps, n = n, K4 = K4)

  additional_term_skewness <-
    0.037 * e_1n(eps = eps, noskewness = FALSE) * lambda3^2 / n

  remainder_term <-
    r2n_inid_skew(eps = eps, n = n, lambda3 = lambda3,
                  K3tilde = K3tilde, K4 = K4, K3 = K3)

  return(main_term + additional_term_skewness + remainder_term)
}

#' Compute the bound of Theorem A.4
#' (Theorem 3.1 under Assumption 2.1: inid case)
#' in the noskewness case (Equation (24))
#' WITHOUT the term involving the integral of f_{S_n}

#' @examples
#' Bound_EE1_cont_inid_noskew_wo_int_fSn(n = 300, eps = 0.1, K4 = 9, K3tilde = 6)
#' Bound_EE1_cont_inid_noskew_wo_int_fSn(n = 300, eps = 0.1, K4 = 12, K3tilde = 6)
#'
#' @noRd
Bound_EE1_cont_inid_noskew_wo_int_fSn <- function(n, eps, K4, K3tilde)
{
  main_term <-
    Bound_EE1_cont_common_part_noskewness(eps = eps, n = n, K4 = K4)

  remainder_term <-
    r2n_inid_noskew(eps = eps, n = n, K3tilde = K3tilde, K4 = K4)

  return(main_term + remainder_term)
}

#' Compute the bound of Theorem A.5
#' (Theorem 3.1 under Assumption 2.2: iid case)
#' in the skewness case (Equation (29))
#' WITHOUT the term involving the integral of f_{S_n}
#'
#' @examples
#' Bound_EE1_cont_iid_skew_wo_int_fSn(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, lambda3 = 1, K3 = 4)
#' Bound_EE1_cont_iid_skew_wo_int_fSn(n = 500, eps = 0.1, K4 = 9, K3tilde = 6, lambda3 = 1, K3 = 4)
#'
#' @noRd
#'
Bound_EE1_cont_iid_skew_wo_int_fSn <- function(n, eps, K4, K3, lambda3, K3tilde)
{
  main_term <-
    Bound_EE1_cont_common_part_noskewness(eps = eps, n = n, K4 = K4)

  additional_term_skewness <-
    0.037 * e_3(eps = eps) * lambda3^2 / n

  remainder_term <-
    r2n_iid_skew(eps = eps, n = n, lambda3 = lambda3,
                 K3tilde = K3tilde, K4 = K4, K3 = K3)

  return(main_term + additional_term_skewness + remainder_term)
}

#' Compute the bound of Theorem A.5
#' (Theorem 3.1 under Assumption 2.2: iid case)
#' in the noskewness case (Equation (30))
#' WITHOUT the term involving the integral of f_{S_n}
#'
#' @examples
#' Bound_EE1_cont_iid_noskew(n = 300, eps = 0.1, K4 = 9, K3tilde = 6, K3 = 4)
#' Bound_EE1_cont_iid_noskew(n = 300, eps = 0.02, K4 = 9, K3tilde = 6, K3 = 4)
#'
#' @noRd
#'
Bound_EE1_cont_iid_noskew_wo_int_fSn <- function(n, eps, K4, K3, K3tilde)
{
  main_term <-
    Bound_EE1_cont_common_part_noskewness(eps = eps, n = n, K4 = K4)

  remainder_term <-
    r2n_iid_noskew(eps = eps, n = n, K3tilde = K3tilde, K4 = K4, K3 = K3)

  return(main_term + remainder_term)
}

#------------------------------------------------------------------------------#
#      Main term                                                           #####
#------------------------------------------------------------------------------#

#' Main term of Theorem A.4 and A.5, excluding skewness
#'
#' @examples
#' Bound_EE1_cont_common_part_noskewness(n = 300, eps = 0.1, K4 = 9)
#' Bound_EE1_cont_common_part_noskewness(n = 1000, eps = 0.1, K4 = 9)
#'
#' @noRd
#'
Bound_EE1_cont_common_part_noskewness <- function(n, eps, K4){
  return(
    0.327 * K4 * (1/12 + 1 / (4 * (1 - 3*eps)^2)) / n
  )
}

#------------------------------------------------------------------------------#
#      Remainder terms                                                     #####
#------------------------------------------------------------------------------#

#' Compute the remainder r_{2,n}^{inid, skew}(epsilon)
#' defined in Equation (27) of the paper, in section A.4
#'
#' @examples
#' r2n_inid_skew(eps = 0.1, n = 200, lambda3 = 2, K3tilde = 6, K3 = 5, K4 = 9)
#' r2n_inid_skew(eps = 0.1, n = 500, lambda3 = 2, K3tilde = 6, K3 = 5, K4 = 12)
#'
#' @noRd
#'
r2n_inid_skew <- function(n, eps, K4, K3, lambda3, K3tilde)
{
  value_Rn_inid_integrated <- Rn_inid_integrated(
    eps = eps, n = n, K4 = K4, lambda3 = lambda3, noskewness = FALSE)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r2n(
    eps = eps, p = 3, n = n, K4 = K4, K3tilde = K3tilde)

  upper_end_Gamma <- 16 *pi^3 * n^2 / K3tilde^4
  lower_end_Gamma <- min( sqrt(2*eps) * (n/K4)^(1/4) , upper_end_Gamma )

  return(
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
      0.3334 * K3tilde^4 * abs(lambda3) / (16 * pi^4 * n^(5/2)) +
      14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
      4.3394 * K3tilde^12 * abs(lambda3) / ((2*pi)^12 * n^(13/2)) +
      abs(lambda3) * abs(Upper_Incomplete_gamma(3/2, upper_end_Gamma) -
                           Upper_Incomplete_gamma(3/2, lower_end_Gamma) ) / sqrt(n) +
      value_Rn_inid_integrated +
      Value_cst_bound_modulus_psi() * K3 * value_Delta_curly_brace_part / (6 * pi * sqrt(n)) +
      common_diffGamma_r2n(n = n, K3tilde = K3tilde)
  )
}

#' Compute the remainder r_{2,n}^{inid, noskew}(epsilon)
#' defined in Equation (28) of the paper, in section A.4
#'
#' @examples
#' r2n_inid_noskew(eps = 0.1, n = 200, K3tilde = 6, K4 = 9)
#' r2n_inid_noskew(eps = 0.1, n = 500, K3tilde = 6, K4 = 12)
#'
#' @noRd
#'
r2n_inid_noskew <- function(eps, n, K3tilde, K4)
{
  value_Rn_inid_integrated <- Rn_inid_integrated(
    eps = eps, noskewness = TRUE, n = n, K4 = K4)

  value_Delta_curly_brace_part <- Delta_curly_brace_part_r2n(
    eps = eps, p = 4, n = n, K4 = K4, K3tilde = K3tilde)

  return(
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
      14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
      value_Rn_inid_integrated +
      Value_cst_bound_modulus_psi() * K4 * value_Delta_curly_brace_part / (6 * pi * n) +
      common_diffGamma_r2n(n = n, K3tilde = K3tilde)
  )
}

#' Compute the remainder r_{2,n}^{iid, skew}(epsilon)
#' defined in Equation (33) of the paper, in section A.5
#'
#' @examples
#' r2n_iid_skew(eps = 0.1, n = 200, lambda3 = 2, K3tilde = 6, K3 = 5, K4 = 9)
#' r2n_iid_skew(eps = 0.1, n = 200, lambda3 = 2, K3tilde = 6, K3 = 5, K4 = 12)
#'
#' @noRd
#'
r2n_iid_skew <- function(eps, n, lambda3, K3tilde, K4, K3)
{
  value_Rn_iid_integrated <- Rn_iid_integrated(
    eps = eps, n = n, K4 = K4, lambda3 = lambda3, noskewness = FALSE)

  upper_end_Gamma <- 16 * pi^3 * n^2 / K3tilde^4
  lower_end_Gamma <- min( sqrt(2*eps) * (n/K4)^(1/4) , upper_end_Gamma )

  return(
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
      0.3334 * K3tilde^4 * abs(lambda3) / (16 * pi^4 * n^(5/2)) +
      14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
      4.3394 * K3tilde^12 * abs(lambda3) / ((2*pi)^12 * n^(13/2)) +
      abs(lambda3) * (Upper_Incomplete_gamma(3/2, lower_end_Gamma) -
                        Upper_Incomplete_gamma(3/2, upper_end_Gamma)) / sqrt(n) +
      value_Rn_iid_integrated +
      1.306 * ( e_2n(eps = eps, noskewness = FALSE , n = n, K4 = K4, lambda3 = lambda3)
                - e_3(eps = eps) ) * lambda3^2 / (36 * n) +
      common_diffGamma_r2n(n = n, K3tilde = K3tilde)
  )
}
# TODO: manque le bout en 1.0253 x ??

#' Compute the remainder r_{2,n}^{iid, noskew}(epsilon)
#' defined in Equation (34) of the paper, in section A.5
#'
#' @examples
#' r2n_iid_noskew(eps = 0.1, n = 200, K3tilde = 6, K4 = 9, K3 = 4)
#' r2n_iid_noskew(eps = 0.1, n = 500, K3tilde = 6, K4 = 12, K3 = 7)
#'
#' @noRd
#'
r2n_iid_noskew <- function(eps, n, K3tilde, K4, K3)
{
  value_Rn_iid_integrated <- Rn_iid_integrated(
    eps = eps, n = n, K4 = K4, noskewness = TRUE)

  upper_end_Gamma <- 2^5 * pi^6 * n^4 / K3^8
  lower_end_Gamma <- min( eps * sqrt(n / (16 * K4)) , upper_end_Gamma)

  return(
    1.2533 * K3tilde^4 / (16 * pi^4 * n^2) +
      14.1961 * K3tilde^16 / ((2*pi)^16 * n^8) +
      value_Rn_iid_integrated +
      16 * Value_cst_bound_modulus_psi() * K3 *
      abs( Upper_Incomplete_gamma(2, upper_end_Gamma) -
             Upper_Incomplete_gamma(2, lower_end_Gamma) ) /(3 * pi * n) +
      common_diffGamma_r2n(n = n, K3tilde = K3tilde)
  )
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
    upper_end = Delta * min(2 * eps * sqrt(n/K4) ,
                            2^8 * pi^6 * n^4 / K3tilde^8)
    lower_end = 2^8 * pi^6 * Delta * n^4 / K3tilde^8

    value <- 0.5 * abs(Delta)^(- p/2) *
      abs( Lower_Incomplete_gamma(p / 2, lower_end) -
             Lower_Incomplete_gamma(p / 2, upper_end) )
  }

  return(value)
}

#' @examples
#' common_diffGamma_r2n(n = 300, K3tilde = 6)
#'
#' @noRd
common_diffGamma_r2n <- function(n, K3tilde)
{
  bound_modulus_psi <- Value_cst_bound_modulus_psi()
  t1star <- Value_t1star()

  T <- 16 * pi^4 * n^2 / K3tilde^4

  shortcut <- (1 - 4 * pi * Value_chi1() * t1star)

  J4 <- bound_modulus_psi / pi *
    abs(Upper_Incomplete_gamma(0, min(T^(1/2), T^2) * shortcut / (2 * pi^2) ) -
          Upper_Incomplete_gamma(0, min(t1star^2 * T^(1/2), T^2 / pi^2) * shortcut / 2) )

  J5 <- bound_modulus_psi / pi *
    abs(Upper_Incomplete_gamma(0, min(T^(1/2), T^2) / (2 * pi^2) ) -
          Upper_Incomplete_gamma(0, T^2 /(2 * pi^2) ) )

  return (J4 + J5)
}


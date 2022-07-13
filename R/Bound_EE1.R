

#' Uniform bound on Edgeworth expansion
#'
#' This function computes a non-aymptotically uniform bound on
#' the difference between a cdf and its 1st order Edgeworth expansion
#'
#' @param setup logical vector of size 3 made up of
#' the following components: \itemize{
#'    \item \code{continuity}: if \code{TRUE}, assume that the distribution is continuous.
#'
#'    \item \code{iid}: if \code{TRUE}, assume that the random variables are i.i.d.
#'
#'    \item \code{no_skewness}: if \code{TRUE}, assume that the distribution is unskewed.
#' }
#'
#' @param regularity list of length up to three
#' (only used in the `continuity=TRUE` framework)
#' with the following components:\itemize{
#'
#'    \item `C0` and `p`: only used in the `iid=FALSE` case.
#'    It corresponds to the assumption of a polynomial bound on f_Sn:
#'    \eqn{|f_Sn(u)| <= C_0 * u^(-p)} for every \eqn{u > a_n},
#'    where \eqn{a_n := 2t_1^* \pisqrt(n) / K3tilde}.
#'
#'    \item `kappa`: only used in the `iid=TRUE` case.
#'    Corresponds to a bound on the modulus of the characteristic function of
#'    the standardized \eqn{X_n}. More precisely, `kappa` is an upper bound on
#'
#'    kappa = sup of modulus of f_{X_n / sigma_n}(t)
#'    over all t such that \eqn{|t| >= 2 t_1^* \pi / K3tilde}
#' }
#'
#' @examples
#' setup = list(continuity = TRUE, iid = FALSE, no_skewness = TRUE)
#' regularity = list(C0 = 1, p = 2, kappa = 0.99)
#'
#' computedBound <- Bound_EE1(
#'   setup = setup, n = 150, K4 = 9,
#'   regularity = regularity, eps = 0.1 )
#'
#' print(computedBound)
#'
#' @export
#'
Bound_EE1 <- function(
  setup = list(continuity = FALSE, iid = FALSE, no_skewness = FALSE),
  n, K4, K3 = NULL, lambda3 = NULL, K3tilde = NULL,
  regularity = list(C0 = 1, p = 2, kappa = 0.99),
  eps = 0.1)
{

  if ( !any(sapply(X = setup, FUN = is.logical)) && length(setup) != 3){
    stop("'setup' should be a logical vector of size 3.")
  }

  env <- environment()
  BEE.updateBounds(env)
  continuity = setup$continuity
  iid = setup$iid
  no_skewness = setup$no_skewness

  if (!continuity & !iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_inid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (!continuity & !iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_inid_noskew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (!continuity & iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_iid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (!continuity & iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_iid_noskew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & !iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_inid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & !iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_inid_noskew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_iid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_iid_noskew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  }

  return(ub_DeltanE)
}


#------------------------------------------------------------------------------#
##### Nocontinuity, inid - Theorem A.2 (Thm 2.1 under Assumption 2.1) #####
#------------------------------------------------------------------------------#

r1n_inid_skew <- function(eps, n, lambda3, K3tilde, K4)
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

Bound_EE1_nocont_common_part_noskewness <- function(eps, n, K4, K3tilde)
{
  return(
    0.1995 * K3tilde / sqrt(n) +
      (0.031 * K3tilde^2 + 0.327 * K4 * (1/12 + 1 / (4 * (1 - 3*eps)^2))) / n
  )
}



Bound_EE1_nocont_inid_main_part <- function(
    eps, noskewness, n, K3tilde, K4, lambda3 = NULL)
{

  if ( TRUE ) {

  } else {

    part_additional_skewness <-
      (0.054 * abs(lambda3) * K3tilde +
         0.037 * e_1n(eps = eps, noskewness = noskewness) * lambda3^2) / n

    value <- part_noskewness + part_additional_skewness
  }

  return(value)
}

Bound_EE1_nocont_inid_skew <- function(n, eps, K4, K3tilde, lambda3)
{
  dominant_term <-
    Bound_EE1_nocont_common_part_noskewness(eps = eps, n = n, K4 = K4, K3tilde = K3tilde) +
    (0.054 * abs(lambda3) * K3tilde +
       0.037 * e_1n(eps = eps, noskewness = FALSE) * lambda3^2) / n

  remainder_term <-
    r1n_inid_skew(eps = eps, n = n, lambda3 = lambda3, K3tilde = K3tilde, K4 = K4)

  return(dominant_term + remainder_term)
}

Bound_EE1_nocont_inid_noskew <- function(n, eps, K4, K3tilde)
{
  dominant_term <-
    Bound_EE1_nocont_common_part_noskewness(eps = eps, n = n, K4 = K4, K3tilde = K3tilde)

  remainder_term <- r1n_inid_noskew(eps = eps, n = n, K4 = K4, K3tilde = K3tilde)

  return(dominant_term + remainder_term)
}

#------------------------------------------------------------------------------#
#####  #####
#------------------------------------------------------------------------------#





Bound_EE1_nocont_inid_skew <- function(n, eps, K4, K3, lambda3, K3tilde)
{

  Rn_inid_integrated_of_eps <- Rn_inid_integrated(
    eps = eps, n = n, K4 = K4, lambda3 = lambda3, noskewness = FALSE)




  stop("replace the p by its value")
  stop("fix the \\Gamma and \\gamma")

  min = min( sqrt(2 * eps) * (n/K4)^(1/4) , 2 * sqrt(n) / K3tilde )

  rninidskew = (14.1961 + 67.0415) * K3tilde^4 / (16 * pi^4 * n^2) +
    4.3394 * abs(lambda3) * K3tilde^3 / (8 * pi^3 * n^2) +
    R1int +
    abs(lambda3) *
    ( Gamma( 3/2 , min )
      - Gamma( 3/2 , 2 * sqrt(n) / K3tilde) )  / sqrt(n) +

    1.0253 * K3 / (6 * pi * sqrt(n)) *
    ( 0.5 * abs(Delta)^(-1-3/2) *
      abs( gamma(3/2, 4 * Delta * n / K3tilde^2)
           - gamma(3/2, Delta * min^2 ) ) )

  result = 1.2533 * K3tilde / (2 * pi * sqrt(n)) +
    0.3334 * abs(lambda3) * K3tilde / (2 * pi * n) +
    1.2187 * K3tilde^2 / (4 * pi * n) +
    (0.327 * K4 / n) * (1/12 + 1 / (4 * (1 - 3 * eps)^2 ) ) +
    (1.306 * e_1(eps = eps, noskewness = FALSE) * lambda3^2) / (36 * n) +
    rninidskew

  return (result)
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

    value <- 0.5 * abs(Delta)^(-1 - p/2) * abs(
      Lower_Incomplete_gamma(p / 2, 4 * Delta * n / K3tilde^2) -
        Lower_Incomplete_gamma(p / 2,
                               2 * Delta *
                                 min(eps * sqrt(n/K4), 2 * n / K3tilde^2)))
  }

  return(value)
}



#------------------------------------------------------------------------------#
##### Nocontinuity, iid - Theorem A.3 (Thm 2.1 under Assumption 2.2) #####
#------------------------------------------------------------------------------#

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


Bound_EE1_nocont_iid_skew <- function(n, eps, K4, K3tilde, lambda3)
{
  dominant_term <-
    Bound_EE1_nocont_common_part_noskewness(eps = eps, n = n, K4 = K4, K3tilde = K3tilde) +
    (0.054 * abs(lambda3) * K3tilde +
       0.037 * e3(eps) * lambda3^2) / n

  remainder_term <-
    r1n_iid_skew(eps = eps, n = n, K3tilde = K3tilde, K4 = K4, lambda3 = lambda3)

  return(dominant_term + remainder_term)
}


Bound_EE1_nocont_iid_skew <- function(n, eps, K4, K3tilde, lambda3)
{
  dominant_term <-
    Bound_EE1_nocont_common_part_noskewness(eps = eps, n = n, K4 = K4, K3tilde = K3tilde)

  remainder_term <-
    r1n_iid_noskew(eps = eps, n = n, K3tilde = K3tilde, K4 = K4)

  return(dominant_term + remainder_term)
}

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






#------------------------------------------------------------------------------#
##### #####
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
##### #####
#------------------------------------------------------------------------------#


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



#' Update missing moments based on upper bounds
#'
#' Indeed, K3, lambda3, K3tilde optional since an upper bound on K4 is enough to
#' upper bound them.
#'
#' @param env an environment including a vector named K4.
#'
#' @return NULL.
#' But after running this function,
#' K3, lambda3, and K3tilde are set if they were NULL before.
#'
#' @noRd
#'
BEE.updateBounds <- function(env)
{
  # Bound (by K4) on K3 if its value (or bound on it) is no provided
  if (is.null(env$K3)){
    env$K3 <- env$K4^(3/4)
  }

  # Bound (by K3) on lambda3 if its value (or bound on it) is no provided
  # In fact, a bound on abs(lambda3) and only abs(lambda3) or lambda3^2
  # is involved, thus fine.
  if (is.null(env$lambda3)){
    env$lambda3 <- Value_cst_bound_lambda3_by_K3() * env$K3
  }

  # Bound on K3tilde (by K3) if its value (or bound on it) is no provided
  if (is.null(env$K3tilde)){
    if (env$setup$iid){
      env$K3tilde <- env$K3 + 1
    } else {
      env$K3tilde <- 2 * env$K3
    }
  }

  return (NULL)
}

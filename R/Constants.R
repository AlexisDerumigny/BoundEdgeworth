
# These are the functions to call recurrent numerical constants

#' @return the value of the constant \eqn{\chi_1}
#' See paragraph 'Additional notation' at the end of section 1 of the paper
#'
#' @noRd
#'
Value_chi1 <- function(){return(0.09916191)}

#' @return the value of the constant \eqn{t_1^*}
#' See paragraph 'Additional notation' at the end of section 1 of the paper
#'
#' @noRd
#'
Value_t1star <- function(){return(0.6359658)}

#' @return the value of the numerical constant that appears
#' to upper bound the modulus of the function \eqn{\Psi}
#' See Equation (9), section A.1 of the paper:
#' \eqn{|Psi(t)| \leq 1.0253 / (2 * pi * |t|)}
#' N.B.: the denominator 2 pi is not included as it sometimes cancel out;
#' it only returns the numerator: 1.0253
#'
#' @noRd
#'
Value_cst_bound_modulus_psi <- function(){return(1.0253)}

#' @return the value of the numerical constant that appears
#' to upper bound the absolute value of lambda3 by K3 (Pinelis, 2011, Theorem 1)
#' \eqn{\abs{\lambda_3} \leq 0.621 K_3}
#' See for instance in Example 2.2 of the paper
#'
#' @noRd
#'
Value_cst_bound_abs_lambda3_by_K3 <- function(){return(0.621)}

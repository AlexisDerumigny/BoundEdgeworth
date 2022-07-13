


Value_chi1 <- function(){
  return(0.09916191)
}

Value_t1star <- function(){
  return(0.6359658)
}

# Majoration
# |Psi(t)| <= 1.0253 / (2 * pi * |t|)
# La constante numerique sans le 2 * pi qui parfois se simplifie.
Value_cst_bound_modulus_psi <- function(){
  return(1.0253)
}

# Numerical constant to upper bound lambda3 by K3 (Pinelis, 2011, Theorem 1)
Value_cst_bound_lambda3_by_K3 <- function(){
  return(0.621)
}

Value_chi1 <- function(){
  return(0.09916191)
}


####


#' Wrap-up for the standard (complete) gamma function
#' \Gamma with the notation of the paper
#' Already implemented as such in R.
#'
#' @noRd
#'
Standard_gamma <- function(x){
  gamma(x)
}

#' (Upper) Incomplete gamma function
#' \Gamma(a,x) := \int_x^{+\infty} u^{a-1} e^{-u} du
#' Can be computed numerically using the package \texttt{expint}~\citep{goulet2016expint} in R.
#' Reference: https://search.r-project.org/CRAN/refmans/expint/html/gammainc.html
#' Since we use it only for non-negative value of a, we use
#' Γ(a,x)=Γ(a)(1−P(a,x)),
#' where Γ(a) is the function implemented by R's gamma() and P(a, x)
#' is the cumulative distribution function of the gamma distribution
#' (with scale equal to one) implemented by R's pgamma().
#' By doing so, we avoid dependence on the expint package.
#'
#' @noRd
#'
Upper_Incomplete_gamma <- function(a, x){
  gamma(a) * pgamma(x, a, rate = 1, lower.tail = FALSE)
}
# TODO: erreur, formulation precedente que dans a > 0
# sinon, utiliser expint::gammainc(), ok pour a = 0
# expint::gammainc(a = 0, x = 10)

#' Lower Incomplete Gamma function
#' \gamma(a, x) := \int_0^x |v|^{a-1} \exp(-v) dv.
#'
#' Vérifié par Alexis et Yannick le 11 juillet
#' Bien de 0 à x avec x qui peut être négatif
#' Bien la valeur absolue dans l'intégrale donc quand x est négatif, ce n'est
#' pas une Gamma, c'est une autre fonction
#'
#'
#' @noRd
#'
Lower_Incomplete_gamma <- function(a, x){
  gamma(a) * pgamma(x, a, rate = 1, lower.tail = TRUE)
}

#' TODO: attention pas a proprement parler une gamma dans le cas ou x est
#' negatif (car valeur absolue), donc a integrer numeriquement probablement
#' sauf si existe un package?
#' Intervient pour Delta negatif (tres rare donc ok de toute maniere)

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

  return(A1n + A2n + A3n + A4n + A5n + A6n + A7n)
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

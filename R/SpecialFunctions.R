
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
  gamma(a) * stats::pgamma(x, a, rate = 1, lower.tail = FALSE)
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
  gamma(a) * stats::pgamma(x, a, rate = 1, lower.tail = TRUE)
}

#' TODO: attention pas a proprement parler une gamma dans le cas ou x est
#' negatif (car valeur absolue), donc a integrer numeriquement probablement
#' sauf si existe un package?
#' Intervient pour Delta negatif (tres rare donc ok de toute maniere)



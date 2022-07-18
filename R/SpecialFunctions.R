
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
  if (a > 0) {
    return (gamma(a) * stats::pgamma(x, a, rate = 1, lower.tail = FALSE))

  } else {
    return (expint::gammainc(a, x))
  }

}

#' Lower Incomplete Gamma function
#' \gamma(a, x) := \int_0^x |v|^{a-1} \exp(-v) dv.
#'
#'
#' @noRd
#'
Lower_Incomplete_gamma <- function(a, x){
  if (x >= 0) {
    return ( gamma(a) * stats::pgamma(x, a, rate = 1, lower.tail = TRUE) )
  } else {
    stop("not implemented yet")
  }
}

#' TODO: attention pas a proprement parler une gamma dans le cas ou x est
#' negatif (car valeur absolue), donc a integrer numeriquement probablement
#' sauf si existe un package?
#' Intervient pour Delta negatif (tres rare donc ok de toute maniere)



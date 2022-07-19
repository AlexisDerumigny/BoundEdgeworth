
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
    return (suppressWarnings( { expint::gammainc(a, x) } ) )
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
    return ( Generalized_Lower_Incomplete_gamma(a = a, x = x) )
  }
}

#' Generalized_Lower_Incomplete_gamma
#'
#' Compute \eqn{\gamma(a,x)} for negative x
#'
#' @noRd
#'
Generalized_Lower_Incomplete_gamma <- function(a, x){

  f_integrand <- Create_Integrand(a)

  res_integrate <- stats::integrate(f = f_integrand,
                                    lower = 0, upper = x)

  return(res_integrate$value)
}

#' Factory function to create the integrand used in
#' Generalized_Lower_Incomplete_gamma
#'
#' @noRd
#'
Create_Integrand <- function(a){
  function(u){
    abs(u)^(a-1) * exp(u)
  }
}


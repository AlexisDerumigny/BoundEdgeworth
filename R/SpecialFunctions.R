
# These are special functions used in our bounds
# See paragraph 'Additional notation' at the end of section 1 of the paper

#' Wrap-up for the standard (complete) gamma function
#' denoted \eqn{\Gamma} in the paper.
#' Already implemented as such in R but to avoid confusion with the lower
#' and upper incomplete Gamma functions.
#'
#' @noRd
#'
Standard_gamma <- function(x){
  return( gamma(x) )
}

#' Upper Incomplete gamma function
#' denoted \eqn{\Gamma(a,x) := \int_x^{+\infty} u^{a-1} e^{-u} du} in the paper.
#' Can be computed numerically using the package \texttt{expint}~\citep{goulet2016expint} in R.
#' Reference: https://search.r-project.org/CRAN/refmans/expint/html/gammainc.html
#'
#' @noRd
#'
Upper_Incomplete_gamma <- function(a, x){
  return( suppressWarnings( { expint::gammainc(a = a, x = x) } ) )
}

#' Lower Incomplete (generalized) Gamma function
#' \gamma(a, x) := \int_0^x |v|^{a-1} \exp(-v) dv.
#' generalized since if x is negative, it is not the usual lower incomplete
#' gamma function due to the absolute value in the integrand.
#'
#' @examples
#' Lower_Incomplete_gamma(a = 2, x = 3)
#' Lower_Incomplete_gamma(a = 2, x = -0.5)
#' # Lower_Incomplete_gamma(a = 0, x = -0.5) # defined for a > 0 only
#' Lower_Incomplete_gamma(a = 2, x = 0) # ok equal to 0
#'
#' @noRd
#'
Lower_Incomplete_gamma <- function(a, x){
  if (x >= 0) {
    return( gamma(a) * stats::pgamma(x, a, rate = 1, lower.tail = TRUE) )
  } else {
    return( Generalized_Lower_Incomplete_gamma(a = a, x = x) )
  }
}

#' Compute \eqn{\gamma(a,x)} for negative x
#' In this case, the function is not the usual lower incomplete gamma function;
#' we compute it by numerical integration using \code{stats::integrate}.
#'
#' @noRd
#'
Generalized_Lower_Incomplete_gamma <- function(a, x){

  f_integrand <- Create_Integrand(a = a)

  res_integrate <- stats::integrate(f = f_integrand, lower = 0, upper = x)

  return( res_integrate$value )
}

#' Factory function to create the integrand used in
#' Generalized_Lower_Incomplete_gamma()
#'
#' @noRd
#'
Create_Integrand <- function(a){
  return( function(u){abs(u)^(a - 1) * exp(-u)} )
}

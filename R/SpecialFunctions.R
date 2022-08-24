
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
Upper_incomplete_gamma <- function(a, x){
  return( suppressWarnings( { expint::gammainc(a = a, x = x) } ) )
}

#' (extended) Lower Incomplete Gamma function
#' \gamma(a, x) := \int_0^x |v|^{a-1} \exp(-v) dv.
#' extended since if x is negative, it is not the usual lower incomplete
#' gamma function due to the absolute value in the integrand.
#'
#' @examples
#' Lower_incomplete_gamma(a = 2, x = 3)
#' Lower_incomplete_gamma(a = 2, x = -0.5)
#' # Lower_incomplete_gamma(a = 0, x = -0.5) # defined for a > 0 only
#' Lower_incomplete_gamma(a = 2, x = 0) # ok equal to 0
#'
#' @noRd
#'
Lower_incomplete_gamma <- function(a, x){
  if (x >= 0) {
    return( gamma(a) * stats::pgamma(x, shape = a, rate = 1, lower.tail = TRUE) )
  } else {
    return( Lower_incomplete_gamma_for_negative_x(a = a, x = x) )
  }
}

#' Compute \eqn{\gamma(a,x)} for negative x
#' In this case, the function is not the usual lower incomplete gamma function;
#' we compute it by numerical integration using \code{stats::integrate}.
#'
#' @noRd
#'
Lower_incomplete_gamma_for_negative_x <- function(a, x)
{
  # Function to be integrated
  f_integrand <- function(u){abs(u)^(a - 1) * exp(-u)}

  res_integrate <- stats::integrate(f = f_integrand, lower = 0, upper = x)

  return( res_integrate$value )
}

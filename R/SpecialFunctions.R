
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
#' stats::pgamma(q = 3, shape = 2, rate = 1, lower.tail = TRUE)
#'
#' Lower_incomplete_gamma(a = 2, x = -0.5)
#' Lower_incomplete_gamma_for_negative_x(a = 2, x = -0.5)
#'
#' # Lower_incomplete_gamma(a = 0, x = -0.5) # defined for a > 0 only
#'
#' Lower_incomplete_gamma(a = 2, x = c(-0.5, -0.2, 3, 8))
#'
#' @noRd
#'
Lower_incomplete_gamma <- function(a, x){
  result <- rep(0, length(x))

  which_positive <- which(x > 0)
  if (length(which_positive > 0)){
    result[which_positive] <- gamma(a) * stats::pgamma(x[which_positive], shape = a,
                                                       rate = 1, lower.tail = TRUE)
  }

  which_negative <- which(x < 0)
  if (length(which_negative > 0)){
    result[which_negative] <- Lower_incomplete_gamma_for_negative_x(a = a, x = x[which_negative])
  }

  return (result)
}

#' Compute \eqn{\gamma(a,x)} for negative x
#' In this case, the function is not the usual lower incomplete gamma function;
#' we compute it by numerical integration using \code{stats::integrate}.
#'
#' @examples
#' # We check that it is the same
#' Lower_incomplete_gamma_for_negative_x(a = 2, x = -0.5)
#'
#' f_integrand_here <- function(u){abs(u)^(2 - 1) * exp(-u)}
#' stats::integrate(f = f_integrand_here, lower = 0, upper = -0.5)
#'
#' @noRd
#'
Lower_incomplete_gamma_for_negative_x <- function(a, x)
{
  # Function to be integrated
  f_integrand <- function(u){abs(u)^(a - 1) * exp(-u)}

  res_integrate <- unlist(lapply(X = x, FUN = function (x) {
    (stats::integrate(f = f_integrand, lower = 0, upper = x))$value} ) )

  return( res_integrate )
}

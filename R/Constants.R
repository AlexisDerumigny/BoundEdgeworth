

#' Value_chi1()
#'
#' @return Return the value of the constant \eqn{\chi_1}
#'
#' @noRd
#'
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


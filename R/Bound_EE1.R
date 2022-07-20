

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
#' @param regularity list of length up to 3
#' (only used in the `continuity=TRUE` framework)
#' with the following components:\itemize{
#'
#'    \item `C0` and `p`: only used in the `iid=FALSE` case.
#'    It corresponds to the assumption of a polynomial bound on f_Sn:
#'    \eqn{|f_Sn(u)| <= C_0 * u^(-p)} for every \eqn{u > a_n},
#'    where \eqn{a_n := 2t_1^* \pi \sqrt(n) / K3tilde}.
#'
#'    \item `kappa`: only used in the `iid=TRUE` case.
#'    Corresponds to a bound on the modulus of the characteristic function of
#'    the standardized \eqn{X_n}. More precisely, `kappa` is an upper bound on
#'    kappa = sup of modulus of f_{X_n / sigma_n}(t)
#'    over all t such that \eqn{|t| >= 2 t_1^* \pi / K3tilde}
#' }
#'
#' @examples
#' setup = list(continuity = TRUE, iid = FALSE, no_skewness = TRUE)
#' regularity = list(C0 = 1, p = 2)
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
  n,
  K4 = 9, K3 = NULL, lambda3 = NULL, K3tilde = NULL,
  regularity = list(C0 = 1, p = 2),
  eps = 0.1)
{

  # Check 'setup' argument and define shortcuts
  if ( !all(sapply(X = setup, FUN = is.logical)) || length(setup) != 3) {
    stop("'setup' should be a logical vector of size 3.")
  }
  continuity <- setup$continuity
  iid <- setup$iid
  no_skewness <- setup$no_skewness

  # A bound on K4 needs to be provided.
  # If bounds on lambda3, K3, and K3tilde are not provided,
  # we use upper bounds that can be defined from K4 only.
  env <- environment(); Update_bounds_on_moments(env)

  # No continuity case (moment condition only)
  if (!continuity) {

    if (!iid & !no_skewness) {
      ub_DeltanE = Bound_EE1_nocont_inid_skew (
        n = n, eps = eps, K4 = K4, K3 = K3, K3tilde = K3tilde, lambda3 = lambda3)
    } else if (!iid & no_skewness) {
      ub_DeltanE = Bound_EE1_nocont_inid_noskew (
        n = n, eps = eps, K4 = K4, K3tilde = K3tilde)
    } else if (iid & !no_skewness) {
      ub_DeltanE = Bound_EE1_nocont_iid_skew (
        n = n, eps = eps, K4 = K4, K3 = K3, K3tilde = K3tilde, lambda3 = lambda3)
    } else if (iid & no_skewness) {
      ub_DeltanE = Bound_EE1_nocont_iid_noskew (
        n = n, eps = eps, K4 = K4, K3tilde = K3tilde)
    }

  # Continuity case (additional regularity conditions)
  } else {

    if (!iid & !no_skewness) {
      ub_DeltanE_wo_int_fSn = Bound_EE1_cont_inid_skew_wo_int_fSn (
        n = n, eps = eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
    } else if (!iid & no_skewness) {
      ub_DeltanE_wo_int_fSn = Bound_EE1_cont_inid_noskew_wo_int_fSn (
        n = n, eps = eps, K4 = K4, K3tilde = K3tilde)
    } else if (iid & !no_skewness) {
      ub_DeltanE_wo_int_fSn = Bound_EE1_cont_iid_skew_wo_int_fSn (
        n = n, eps = eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
    } else if (iid & no_skewness) {
      ub_DeltanE_wo_int_fSn = Bound_EE1_cont_iid_noskew_wo_int_fSn (
        n = n, eps = eps, K4 = K4, K3 = K3, K3tilde = K3tilde)
    }

    ub_DeltanE = ub_DeltanE_wo_int_fSn +
      Smoothness_additional_term(n = n, K3tilde = K3tilde,
                                 regularity = regularity, iid = iid)
  }

  return(ub_DeltanE)
}

#' Additional smoothness term
#'
#' It checks if \code{regularity} is well-formated, and if so,
#' it computes the additional smoothness term.
#'
#' @examples
#'
#' Smoothness_additional_term(n = 1000, K3tilde = 6, regularity = list(C0 = 1, p = 2), iid = FALSE)
#' Smoothness_additional_term(n = 1000, K3tilde = 6, regularity = list(kappa = 0.99), iid = TRUE)
#'
#' @noRd
#'
Smoothness_additional_term <- function(n, K3tilde, regularity, iid){

  a_n <- min(2 * Value_t1star() * pi * sqrt(n) / K3tilde,
             16 * pi^3 * n^2 / K3tilde^4 )

  b_n <- 16 * pi^4 * n^2 / K3tilde^4

  success = TRUE

  switch (as.character(length(regularity)),
    "1" = {
      if ((names(regularity) == "kappa") && (iid)) {
        result = Value_cst_bound_modulus_psi() / pi *
          regularity$kappa^n * log(b_n / a_n)
      } else {
        success = FALSE
      }
    },

    "2" = {
      if (identical(names(regularity), c("C0", "p")) ||
          identical(names(regularity), c("p", "C0"))) {
        result = Value_cst_bound_modulus_psi() / pi *
          regularity$C0 * (a_n^(- regularity$p) -  b_n^(- regularity$p))
      } else {
        success = FALSE
      }
    },

    { success = FALSE }
  )

  if (!success) {
    stop("'regularity' should be either a list with C0 and p, or, in the iid case only, a list with only kappa")
  }

  return (result)
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
Update_bounds_on_moments <- function(env) {

  # Bound (by K4) on K3 if its value (or bound on it) is no provided
  if (is.null(env$K3)){
    env$K3 <- env$K4^(3/4)
  }

  # Bound (by K3) on lambda3 if its value (or bound on it) is no provided
  # In fact, a bound on abs(lambda3) and only abs(lambda3) or lambda3^2
  # is involved, thus fine.
  if (is.null(env$lambda3)){
    env$lambda3 <- Value_cst_bound_abs_lambda3_by_K3() * env$K3
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



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
#' @param regularity list of length up to three
#' (only used in the `continuity=TRUE` framework)
#' with the following components:\itemize{
#'
#'    \item `C0` and `p`: only used in the `iid=FALSE` case.
#'    It corresponds to the assumption of a polynomial bound on f_Sn:
#'    \eqn{|f_Sn(u)| <= C_0 * u^(-p)} for every \eqn{u > a_n},
#'    where \eqn{a_n := 2t_1^* \pisqrt(n) / K3tilde}.
#'
#'    \item `kappa`: only used in the `iid=TRUE` case.
#'    Corresponds to a bound on the modulus of the characteristic function of
#'    the standardized \eqn{X_n}. More precisely, `kappa` is an upper bound on
#'
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
  n, K4, K3 = NULL, lambda3 = NULL, K3tilde = NULL,
  regularity = list(C0 = 1, p = 2, kappa = 0.99),
  eps = 0.1)
{

  if ( !any(sapply(X = setup, FUN = is.logical)) && length(setup) != 3){
    stop("'setup' should be a logical vector of size 3.")
  }

  env <- environment()
  BEE.updateBounds(env)
  continuity = setup$continuity
  iid = setup$iid
  no_skewness = setup$no_skewness

  if (!continuity & !iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_inid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (!continuity & !iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_inid_noskew (
      n = n, eps, K4 = K4, K3tilde = K3tilde)
  } else if (!continuity & iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_iid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (!continuity & iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_iid_noskew (
      n = n, eps, K4 = K4, K3tilde = K3tilde)
  } else if (continuity & !iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_inid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & !iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_inid_noskew (
      n = n, eps, K4 = K4, K3 = K3, K3tilde = K3tilde)
  } else if (continuity & iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_iid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_iid_noskew (
      n = n, eps, K4 = K4, K3 = K3, K3tilde = K3tilde)
  }

  if (continuity){
    ub_DeltanE = ub_DeltanE + smoothness_addit_term(n = n, K3tilde = K3tilde,
                                                    regularity = regularity)
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
#' smoothness_addit_term(n = 1000, K3tilde = 6, regularity = list(C0 = 1, p = 2))
#' smoothness_addit_term(n = 1000, K3tilde = 6, regularity = list(kappa = 0.99))
#'
#' @noRd
smoothness_addit_term <- function(
    n, K3tilde,
    regularity = list(C0 = 1, p = 2, kappa = 0.99))
{
  an = min(2 * Value_t1star() * pi * sqrt(n) / K3tilde,
           16 * pi^3 * n^2 / K3tilde^4 )

  bn = 16 * pi^4 * n^2 / K3tilde^4


  success = TRUE

  switch (as.character(length(regularity)),
    "1" = {
      if (names(regularity) == "kappa") {
        result = Value_cst_bound_modulus_psi() * regularity$kappa^n * (bn / an) / pi
      } else {
        success = FALSE
      }
    },

    "2" = {
      if (identical(names(regularity), c("C0", "p")) |
          identical(names(regularity), c("p", "C0")) ) {
        result = Value_cst_bound_modulus_psi() *
          ( regularity$C0 * an^(- regularity$p) - regularity$C0 * bn^(- regularity$p) ) / pi
      } else {
        success = FALSE
      }
    },

    { success = FALSE }
  )

  if (!success){
    stop("'regularity' should be either a list with C0 and p, or a list with only kappa.")
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
BEE.updateBounds <- function(env)
{
  # Bound (by K4) on K3 if its value (or bound on it) is no provided
  if (is.null(env$K3)){
    env$K3 <- env$K4^(3/4)
  }

  # Bound (by K3) on lambda3 if its value (or bound on it) is no provided
  # In fact, a bound on abs(lambda3) and only abs(lambda3) or lambda3^2
  # is involved, thus fine.
  if (is.null(env$lambda3)){
    env$lambda3 <- Value_cst_bound_lambda3_by_K3() * env$K3
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



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
#' regularity = list(C0 = 1, p = 2, kappa = 0.99)
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
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (!continuity & iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_iid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (!continuity & iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_nocont_iid_noskew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & !iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_inid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & !iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_inid_noskew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & iid & !no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_iid_skew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  } else if (continuity & iid & no_skewness)
  {
    ub_DeltanE = Bound_EE1_cont_iid_noskew (
      n = n, eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde)
  }

  return(ub_DeltanE)
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

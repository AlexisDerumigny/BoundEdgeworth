

#' Uniform bound on Edgeworth expansion
#'
#' This function computes a non-aymptotically uniform bound on
#' the difference between the cdf of a normalized sum of random varialbles
#' and its 1st order Edgeworth expansion.
#' It returns a valid value \mjseqn{\delta_n} such that
#' \mjtdeqn{\sup_{x \in R}
#' \left| \textrm{Prob}(S_n \leq x) - \Phi(x)
#' - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) \right|
#' \leq \delta_n,}{
#' \sup_{x \in \mathbb{R}}
#' \left| \textrm{Prob}(S_n \leq x) - \Phi(x)
#' - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) \right|
#' \leq \delta_n,}{
#' \sup_{x \in R} | Prob(S_n \leq x) - \Phi(x)
#' - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) |
#' \leq \delta_n,}
#' where \mjseqn{X_1, \dots, X_n} be \mjseqn{n} independent centered variables,
#' and \mjseqn{S_n} be their normalized sum, in the sense that
#' \mjseqn{S_n := \sum_{i=1}^n X_i / \textrm{sd}(\sum_{i=1}^n X_i)}.
#' Here \mjseqn{\lambda_{3,n}} denotes the average skewness of
#' the variables \mjseqn{X_1, \dots, X_n}.
#'
#' \loadmathjax
#'
#' Note that the variables \mjseqn{X_1, \dots, X_n} must be independent
#' but may have different distributions (if \code{setup$iid = FALSE}).
#'
#'
#' @import mathjaxr
#'
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
#' (only used in the \code{continuity=TRUE} framework)
#' with the following components:\itemize{
#'
#'    \item \code{C0} and \code{p}: only used in the \code{iid=FALSE} case.
#'    It corresponds to the assumption of a polynomial bound on \mjseqn{f_{S_n}}:
#'    \mjseqn{|f_{S_n}(u)| \leq C_0 \times u^{-p}} for every \mjseqn{u > a_n},
#'    where \mjseqn{a_n := 2 t_1^* \pi \sqrt{n} / K3tilde}.
#'
#'    \item \code{kappa}: only used in the \code{iid=TRUE} case.
#'    Corresponds to a bound on the modulus of the characteristic function of
#'    the standardized \mjseqn{X_n}. More precisely, \code{kappa} is an upper bound on
#'    \mjseqn{kappa :=} sup of modulus of \mjseqn{f_{X_n / \sigma_n}(t)}
#'    over all \mjseqn{t} such that \mjseqn{|t| \geq 2 t_1^* \pi / K3tilde}.
#' }
#'
#' @param n sample size ( = number of random variables that appear in the sum).
#'
#' @param K4 bound on the 4th normalized moment of the random variables.
#' We advise to use K4 = 9 as a general case which covers most ``usual'' distributions.
#'
#' @param K3 bound on the 3rd normalized moment.
#' If not given, an upper bound on \code{K3} will be derived from the value of \code{K4}.
#'
#' @param lambda3 (average) skewness of the variables.
#' If not given, an upper bound on \mjseqn{abs(lambda3)}
#' will be derived from the value of \code{K4}.
#'
#' @param K3tilde value of
#' \mjtdeqn{K_{3,n} + \frac{1}{n}\sum_{i=1}^n
#' E|X_i| \sigma_{X_i}^2 / \overline{B}_n^3}{
#' K_{3,n} + \frac{1}{n}\sum_{i=1}^n
#' \mathbb{E}|X_i| \sigma_{X_i}^2 / \overline{B}_n^3}{
#' K_{3,n} + \frac{1}{n}\sum_{i=1}^n E|X_i| \sigma_{X_i}^2 / \overline{B}_n^3}
#' where \mjseqn{\overline{B}_n := \sqrt{(1/n) \sum_{i=1}^n E[X_i^2]}}.
#' If not given, an upper bound on \code{K3tilde} will be derived
#' from the value of \code{K4}.
#'
#' @param eps a value between 0 and 1/3 on which several terms depends.
#' Any value of \code{eps} will give a valid upper bound but some may give
#' tighter results than others.
#'
#' @param verbose if 0 the function is silent (no printing).
#' Higher values gives more precise information about the computation.
#'
#' @return A vector of the same size as \code{n} with values \mjseqn{\delta_n}
#' such that
#' \mjtdeqn{\sup_{x \in R}
#' \left| \textrm{Prob}(S_n \leq x) - \Phi(x)
#' - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) \right|
#' \leq \delta_n.}{
#' \sup_{x \in \mathbb{R}}
#' \left| \textrm{Prob}(S_n \leq x) - \Phi(x)
#' - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) \right|
#' \leq \delta_n.}{
#' \sup_{x \in R} | Prob(S_n \leq x) - \Phi(x)
#' - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) |
#' \leq \delta_n.}
#'
#' @references
#' Derumigny A., Girard L., and Guyonvarch Y. (2021).
#' Explicit non-asymptotic bounds for the distance to the first-order Edgeworth expansion,
#' ArXiv preprint \href{https://arxiv.org/abs/2101.05780}{arxiv:2101.05780}.
#'
#' @seealso \code{\link{Bound_BE}()} for a Berry-Esseen bound.
#'
#'
#' @examples
#' setup = list(continuity = TRUE, iid = FALSE, no_skewness = TRUE)
#' regularity = list(C0 = 1, p = 2)
#'
#' computedBound <- Bound_EE1(
#'   setup = setup, n = c(150, 2000), K4 = 9,
#'   regularity = regularity, eps = 0.1 )
#'
#' setup = list(continuity = TRUE, iid = TRUE, no_skewness = TRUE)
#' regularity = list(kappa = 0.99)
#'
#' computedBound2 <- Bound_EE1(
#'   setup = setup, n = c(150, 2000), K4 = 9,
#'   regularity = regularity, eps = 0.1 )
#'
#' setup = list(continuity = FALSE, iid = FALSE, no_skewness = TRUE)
#'
#' computedBound3 <- Bound_EE1(
#'   setup = setup, n = c(150, 2000), K4 = 9, eps = 0.1 )
#'
#' setup = list(continuity = FALSE, iid = TRUE, no_skewness = TRUE)
#'
#' computedBound4 <- Bound_EE1(
#'   setup = setup, n = c(150, 2000), K4 = 9, eps = 0.1 )
#'
#' print(computedBound)
#' print(computedBound2)
#' print(computedBound3)
#' print(computedBound4)
#'
#' @export
#'
Bound_EE1 <- function(
  setup = list(continuity = FALSE, iid = FALSE, no_skewness = FALSE),
  n,
  K4 = 9, K3 = NULL, lambda3 = NULL, K3tilde = NULL,
  regularity = list(C0 = 1, p = 2),
  eps = 0.1,
  verbose = 0)
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
        n = n, eps = eps, K4 = K4, K3 = K3, lambda3 = lambda3, K3tilde = K3tilde, verbose = verbose)
    } else if (iid & no_skewness) {
      ub_DeltanE_wo_int_fSn = Bound_EE1_cont_iid_noskew_wo_int_fSn (
        n = n, eps = eps, K4 = K4, K3 = K3, K3tilde = K3tilde)
    }

    smoothness_additional_term = Smoothness_additional_term(
      n = n, K3tilde = K3tilde, regularity = regularity, iid = iid)

    if (verbose){
      cat("Smoothness additional term:", smoothness_additional_term, "\n")
    }

    ub_DeltanE = ub_DeltanE_wo_int_fSn + smoothness_additional_term
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

  a_n <- pmin(2 * Value_t1star() * pi * sqrt(n) / K3tilde,
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

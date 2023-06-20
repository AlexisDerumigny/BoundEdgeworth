
#' Computation of uniformly valid power and sufficient sample size for the one-sided Gauss test
#'
#' Let \mjseqn{X_1, \dots, X_n} be \mjseqn{n} i.i.d. variables
#' with mean \mjseqn{\mu}, variance \mjseqn{\sigma^2}.
#' Assume that we want to test the hypothesis
#' \mjseqn{H_0: \mu \leq \mu_0} against the alternative \mjseqn{H_1: \mu \leq \mu_0}.
#' For this, we want to use the classical Gauss test, which rejects the null hypothesis
#' if \mjseqn{\sqrt{n}(\bar{X}_n - \mu)} is larger than the quantile of the Gaussian
#' distribution at level \mjseqn{1 - \alpha}.
#' Let \mjseqn{\eta := (\mu - \mu_0) / \sigma} be the effect size,
#' i.e. the distance between the null and the alternative hypotheses,
#' measured in terms of standard deviations.
#' Let \code{beta} be the uniform power of this test:
#' \mjtdeqn{beta = \inf_{H_1} \textrm{Prob}(\textrm{Rejection}),}{
#' beta = \inf_{H_1} \textrm{Prob}(\textrm{Rejection}),}{
#' beta = \inf_{H_1} \textrm{Prob}(\textrm{Rejection}),}
#' where the infimum is taken over all distributions under the alternative hypothesis, i.e.
#' that have mean \mjseqn{\mu = \mu_0 + \eta \sigma}, bounded kurtosis \code{K4},
#' and that satisfy the regularity condition \code{kappa} described below.
#' This means that this power \code{beta} is uniformly valid over
#' a large (infinite-dimensional) class of alternative distributions,
#' much beyond the Gaussian family even though the test is based on the Gaussian quantile.
#' There is a relation between the sample size \code{n}, the effect size \code{eta}
#' and the uniform power \code{beta} of this test.
#' This function takes as an input two of the three quantities \code{n}, \code{eta}, \code{beta}
#' and return the third one.
#'
#' \loadmathjax
#' @import mathjaxr
#'
#' @param n sample size
#'
#' @param eta the effect size \mjseqn{\eta} that
#' characterizes the alternative hypothesis
#'
#' @param beta the power of detecting the effect \code{eta} using the sample size \code{n}
#'
#' @param alpha the level of the test
#'
#' @param K4 the kurtosis of the \mjseqn{X_i}
#'
#' @param kappa Regularity parameter of the distribution of the \mjseqn{X_i}
#' It corresponds to a bound on the modulus of the characteristic function
#' \mjseqn{f_{X_n / \sigma_n}(t)} of the standardized \mjseqn{X_n}.
#' More precisely, \code{kappa} is an upper bound on
#' \mjseqn{kappa :=} sup of modulus of \mjseqn{f_{X_n / \sigma_n}(t)}
#' over all \mjseqn{t} such that \mjseqn{|t| \geq 2 t_1^* \pi / K3tilde}.
#'
#' @return The computed value of either the sufficient sample size \code{n},
#' or the minimum effect size \code{eta}, or the power \code{beta}.
#'
#' @references
#' Derumigny A., Girard L., and Guyonvarch Y. (2021).
#' Explicit non-asymptotic bounds for the distance to the first-order Edgeworth expansion,
#' ArXiv preprint \href{https://arxiv.org/abs/2101.05780}{arxiv:2101.05780}.
#'
#' @examples
#'
#' # Sufficient sample size to detect an effect of 0.5 standard deviation with probability 80%
#' Gauss_test_powerAnalysis(eta = 0.5, beta = 0.8)
#' # We can detect an effect of 0.5 standard deviations with probability 80% for n >= 548
#'
#' # Power of an experiment to detect an effect of 0.5 with a sample size of n = 800
#' Gauss_test_powerAnalysis(eta = 0.5, n = 800)
#' # We can detect an effect of 0.5 standard deviations with probability 85.1% for n = 800
#'
#' # Smallest effect size that can be detected with a probability of 80% for a sample size of n = 800
#' Gauss_test_powerAnalysis(n = 800, beta = 0.8)
#' # We can detect an effect of 0.114 standard deviations with probability 80% for n = 800
#'
#'
#' @export
#'
Gauss_test_powerAnalysis = function(eta = NULL, n = NULL, beta = NULL,
                                    alpha = 0.05, K4 = 9, kappa = 0.99){

  if (is.null(eta) + is.null(n) + is.null(beta) != 1){
    stop("Exactly two of 'eta', 'n', 'beta' should be known ",
         "( = not set to NULL) to be able to find the third one.")
  }

  if (is.null(beta)){
    result = Gauss_test_power(eta = eta, n = n,
                              alpha = alpha, K4 = K4, kappa = kappa)
    names(result) <- "beta"

  } else if (is.null(n)){

    res = stats::uniroot(
      f = function(n){Gauss_test_power(eta = eta, n = n,
                                       alpha = alpha, K4 = K4, kappa = kappa) - beta},
      interval = c(1, 10^10))

    result = ceiling(res$root)
    names(result) <- "n"

  } else {

    res = stats::uniroot(
      f = function(eta){Gauss_test_power(eta = eta, n = n,
                                         alpha = alpha, K4 = K4, kappa = kappa) - beta},
      interval = c(0, 1))

    result = res[["root"]]
    names(result) <- "eta"
  }

  return (result)
}


# compute beta from eta and n
Gauss_test_power = function(eta, n, alpha, K4, kappa){

  delta_n = BoundEdgeworth::Bound_BE(
    setup = list(continuity = TRUE, iid = TRUE, no_skewness = FALSE),
    n = n, K4 = K4, regularity = list(kappa = kappa), eps = 0.1)

  qnormalpha = stats::qnorm(1 - alpha)

  result = 1 - stats::pnorm( qnormalpha - eta * sqrt(n) ) -
    0.621 * K4^(3/4) * (1 - ( qnormalpha - eta * sqrt(n) )^2 ) *
    stats::dnorm( qnormalpha - eta * sqrt(n) ) / (6 * sqrt(n)) - delta_n

  return (result)
}

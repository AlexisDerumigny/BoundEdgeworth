
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BoundEdgeworth

<!-- badges: start -->
<!-- badges: end -->

This package implements the computation of the bounds described in the
article Derumigny, Girard, and Guyonvarch (2023), Explicit
non-asymptotic bounds for the distance to the first-order Edgeworth
expansion, Sankhya A.
[doi:10.1007/s13171-023-00320-y](https://doi.org/10.1007/s13171-023-00320-y)
[arxiv:2101.05780](https://arxiv.org/abs/2101.05780).

## How to install

You can install the release version from the CRAN:

``` r
install.packages("BoundEdgeworth")
```

or the development version from
[GitHub](https://github.com/AlexisDerumigny/BoundEdgeworth):

``` r
# install.packages("remotes")
remotes::install_github("AlexisDerumigny/BoundEdgeworth")
```

## Available bounds

Let $X_1, \dots, X_n$ be $n$ independent centered variables, and $S_n$
be their normalized sum, in the sense that
$$S_n := \sum_{i=1}^n X_i / \text{sd} \Big(\sum_{i=1}^n X_i \Big).$$

The goal of this package is to compute values of $\delta_n > 0$ such
that bounds of the form

$$
\sup_{x \in \mathbb{R}}
\left| \textrm{Prob}(S_n \leq x) - \Phi(x) \right|
\leq \delta_n,
$$

or of the form

$$
\sup_{x \in \mathbb{R}}
\left| \textrm{Prob}(S_n \leq x) - \Phi(x) - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) \right|
\leq \delta_n,
$$

are valid. Here $\lambda_{3,n}$ denotes the average skewness of the
variables $X_1, \dots, X_n$, $\Phi$ denotes the cumulative distribution
function (cdf) of the standard Gaussian distribution, and $\varphi$
denotes its density.

The first type of bounds is returned by the function `Bound_BE()`
(Berry-Esseen-type bound) and the second type (Edgeworth expansion-type
bound) is returned by the function `Bound_EE1()`.

Such bounds are useful because they can help to control uniformly the
distance between the cdf of a normalized sum $\textrm{Prob}(S_n \leq x)$
and its limit $\Phi(x)$ (which is known by the central limit theorem).
The second type of bound is more precise, and give a control of the
uniform distance between $\textrm{Prob}(S_n \leq x)$ and its first-order
Edgeworth expansion, i.e. the limit from the central limit theorem
$\Phi(x)$ plus the next term
$\frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x)$.

Note that these bounds depends on the assumptions made on
$(X_1, \dots, X_n)$ and especially on $K4$, the average kurtosis of the
variables $X_1, \dots, X_n$. In all cases, they need to have finite
fourth moment and to be independent. To get improved bounds, several
additional assumptions can be added:

- the variables $X_1, \dots, X_n$ are identically distributed,
- the skewness (normalized third moment) of $X_1, \dots, X_n$ are all
  $0$.
- the distribution of $X_1, \dots, X_n$ admits a continuous component.

### Example

``` r
setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE)

Bound_EE1(setup = setup, n = 1000, K4 = 9)
#> [1] 0.1626857
```

This shows that

$$
\sup_{x \in \mathbb{R}}
\left| \textrm{Prob}(S_n \leq x) - \Phi(x) - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) \right|
\leq 0.1626857,
$$

as soon as the variables $X_1, \dots, X_{1000}$ are i.i.d. with a
kurtosis smaller than $9$.

Adding one more regularity assumption on the distribution of the $X_i$
helps to achieve a better bound:

``` r
setup = list(continuity = TRUE, iid = TRUE, no_skewness = FALSE)

Bound_EE1(setup = setup, n = 1000, K4 = 9, regularity = list(kappa = 0.99))
#> [1] 0.1214038
```

This shows that

$$
\sup_{x \in \mathbb{R}}
\left| \textrm{Prob}(S_n \leq x) - \Phi(x) - \frac{\lambda_{3,n}}{6\sqrt{n}}(1-x^2) \varphi(x) \right|
\leq 0.1214038,
$$

in this case.

## Applications to testing

This package also includes the function `Gauss_test_powerAnalysis()`,
that computes a uniformly valid power for the Gauss test that is valid
over a large class of non-Gaussian distribution. This uniform validity
is a consequence of the above-mentioned bounds.

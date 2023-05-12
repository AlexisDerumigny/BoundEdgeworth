
* New function `t_test_powerAnalysis` for computation of power or sufficient sample size
  for the one-sided t-test. The power (respectively sufficient sample size)
  are computed in an exact (non-asymptotic) way based on bounds on the Edgeworth expansion.
  This power holds non-asymptotically and uniformly over all alternatives with a given effect size
  (under some regularity conditions: bounded kurtosis, tail constraint `kappa`).

* Adding a verbose option to `Bound_EE1` to follow the behavior of the different terms.
  Currently, only the setup
  `list(continuity = TRUE, iid = TRUE, no_skewness = FALSE)`
  is covered.

* Removing pipes in `Lower_incomplete_gamma_for_negative_x` to be compatible with R (< 4.1.0).


# BoundEdgeworth 0.1.1

* Documentation update.


# BoundEdgeworth 0.1.0

* Initial release.

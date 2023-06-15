test_that("Bound_EE1 works", {

  for (continuity_ in c(TRUE, FALSE)){
    for (iid_ in c(TRUE, FALSE)){
      for (no_skewness_ in c(TRUE, FALSE)){
        setup = list(continuity = continuity_,
                     iid = iid_,
                     no_skewness = no_skewness_)

        regularity = if (iid_) {list(kappa = 0.99)} else {list(C0 = 1, p = 2)}

        # Test of the verbose option - with printing

        expect_output({
        computedBound <- Bound_EE1(
          setup = setup, n = c(2000), K4 = 9,
          regularity = regularity, eps = 0.1, verbose = 1 )
        })

        # Test of the verbose option - without printing

        expect_silent({
          computedBound <- Bound_EE1(
            setup = setup, n = c(2000), K4 = 9,
            regularity = regularity, eps = 0.1, verbose = 0 )
        })

        # Testing that the bound is positive

        expect_gt(computedBound, 0)
      }
    }
  }

})

test_that("suggest_n_iter runs without error", {

      # nullcat -- example with precomputed trace
      set.seed(1234)
      x <- cat_mat(100)
      expect_no_error({
            trace <- trace_cat(x = x, fun = "nullcat", n_iter = 100,
                               n_chains = 1, method = "curvecat")
            suggest_n_iter(trace, tail_frac = 0.3, plot = FALSE)
      })

      # quantize -- example calling trace internally, and plotting
      set.seed(1234)
      x <- matrix(runif(100), 10)
      expect_no_error(
            n_iter <- suggest_n_iter(
                  x = x, n_chains = 2, n_iter = 200, tail_frac = 0.3,
                  fun = "quantize", n_strata = 4, fixed = "stratum",
                  method = "curvecat", plot = TRUE)
      )

})

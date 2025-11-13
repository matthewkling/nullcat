test_that("`quantize_batch` works", {
      set.seed(123)
      m <- matrix(runif(100), nrow = 10)

      # single core
      expect_no_error(quantize_batch(m, n_reps = 10, n_iter = 50,
                                     n_strata = 3, method = "curvecat"))

      # multicore
      expect_no_error(quantize_batch(m, n_reps = 10, n_iter = 50, n_cores = 4,
                                     n_strata = 3, method = "curvecat"))
})

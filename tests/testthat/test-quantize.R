test_that("quantize works", {
      set.seed(123)
      m <- matrix(runif(100), nrow = 25)

      # prep function runs without error
      expect_no_error(p <- quantize_prep(m, n_iter = 5, n_strata = 3, method = "curvecat", priority = "rows"))

      # main function runs without error, with and without prep supplied
      expect_no_error(quantize(m, prep = p))
      expect_no_error(r <- quantize(m, n_iter = 5, n_strata = 3, method = "curvecat", priority = "rows"))
})


test_that("`priority` correctly preserves marginal sums and conditional distributions", {
      set.seed(123)
      m <- matrix(runif(200), nrow = 20)

      r <- quantize(m, n_iter = 100, n_strata = 3, method = "curvecat", priority = "rows")
      expect_equal(rowSums(m), rowSums(r))
      expect_true(any(colSums(m) != colSums(r)))
      expect_equal(apply(m, 1, sort), apply(r, 1, sort))

      r <- quantize(m, n_iter = 100, n_strata = 3, method = "curvecat", priority = "cols")
      expect_equal(colSums(m), colSums(r))
      expect_true(any(rowSums(m) != rowSums(r)))
      expect_equal(apply(m, 2, sort), apply(r, 2, sort))
})


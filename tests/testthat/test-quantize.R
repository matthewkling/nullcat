test_that("quantize works", {
      set.seed(123)
      m <- matrix(runif(100), nrow = 25)

      # prep function runs without error
      expect_no_error(p <- quantize_prep(m, n_iter = 5, n_strata = 3, method = "curvecat"))

      # main function runs without error, with and without prep supplied
      expect_no_error(quantize(m, prep = p))
      expect_no_error(r <- quantize(m, n_iter = 5, n_strata = 3, method = "curvecat", fixed = "row"))

      # various parameter settings work
      expect_no_error(quantize(m, n_iter = 5, n_strata = 5, method = "curvecat", fixed = "cell"))
      expect_no_error(quantize(m, n_strata = 3, method = "curveball", fixed = "col"))
      expect_no_error(quantize(m, n_iter = 5, breaks = c(.33, .66), method = "curvecat", fixed = "cell"))
})


test_that("`fixed` correctly preserves marginal sums and conditional distributions", {
      set.seed(123)
      m <- matrix(runif(200), nrow = 20)

      r <- quantize(m, n_iter = 100, n_strata = 3, method = "curvecat", fixed = "row")
      expect_equal(rowSums(m), rowSums(r))
      expect_true(any(colSums(m) != colSums(r)))
      expect_equal(apply(m, 1, sort), apply(r, 1, sort))

      r <- quantize(m, n_iter = 100, n_strata = 3, method = "curvecat", fixed = "col")
      expect_equal(colSums(m), colSums(r))
      expect_true(any(rowSums(m) != rowSums(r)))
      expect_equal(apply(m, 2, sort), apply(r, 2, sort))
})


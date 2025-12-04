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
      expect_no_error(quantize(m, n_iter = 5, breaks = c(.33, .66), method = "curvecat", fixed = "cell"))
})


test_that("`fixed` parameter correctly preserves marginal sums and conditional distributions", {
      set.seed(123)
      m <- matrix(runif(200), nrow = 20)

      # "row" should preserve row sums and multisets, but not columns
      fr <- quantize(m, n_iter = 100, n_strata = 3, method = "curvecat", fixed = "row")
      expect_equal(rowSums(m), rowSums(fr))
      expect_true(any(colSums(m) != colSums(fr)))
      expect_equal(apply(m, 1, sort), apply(fr, 1, sort))

      # "col" should preserve column sums and multisets, but not rows
      fc <- quantize(m, n_iter = 100, n_strata = 3, method = "curvecat", fixed = "col")
      expect_equal(colSums(m), colSums(fc))
      expect_true(any(rowSums(m) != rowSums(fc)))
      expect_equal(apply(m, 2, sort), apply(fc, 2, sort))

      # "cell" should preserve neither
      ft <- quantize(m, n_iter = 100, n_strata = 3, method = "curvecat", fixed = "cell")
      expect_false(all(colSums(m) == colSums(ft)))
      expect_false(all(rowSums(m) == rowSums(ft)))
})


test_that("fixed = 'cell' preserves one-to-one mapping with index output", {
      set.seed(1)
      x <- matrix(runif(100), nrow = 10)
      prep <- quantize_prep(x, method = "curvecat", fixed = "cell", n_iter = 1000)

      set.seed(1)
      idx <- curvecat(prep$strata, output = "index")

      set.seed(1)
      xr <- quantize(prep = prep)

      expect_equal(xr, matrix(prep$pool[idx], nrow = nrow(x)))
})


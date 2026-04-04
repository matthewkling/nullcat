test_that("wt_row/wt_col validate correctly", {

      set.seed(123)
      x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

      # Must be numeric matrix
      expect_error(nullcat(x, wt_row = "bad"), "numeric matrix")

      # Must be square
      expect_error(nullcat(x, wt_row = matrix(1, 3, 4)), "square matrix")

      # Must match margin dimension
      expect_error(nullcat(x, wt_row = matrix(1, 5, 5)), "matching nrow")
      expect_error(nullcat(x, wt_col = matrix(1, 5, 5)), "matching ncol")

      # No negative values
      W <- matrix(1, 10, 10)
      W[1, 2] <- -1
      expect_error(nullcat(x, wt_row = W), "non-negative")

      # No NAs
      W <- matrix(1, 10, 10)
      W[1, 2] <- NA
      expect_error(nullcat(x, wt_row = W), "NA")

      # Not supported for r0cat/c0cat
      expect_error(nullcat(x, method = "r0cat", wt_row = matrix(1, 10, 10)),
                   "sequential methods")

      # wt_row incompatible with horizontal swaps
      expect_error(nullcat(x, swaps = "horizontal", wt_row = matrix(1, 10, 10)),
                   "cannot be used with swaps = 'horizontal'")

      # wt_col incompatible with vertical swaps
      expect_error(nullcat(x, swaps = "vertical", wt_col = matrix(1, 10, 10)),
                   "cannot be used with swaps = 'vertical'")

      # Both supplied requires alternating
      expect_error(nullcat(x, swaps = "vertical",
                           wt_row = matrix(1, 10, 10),
                           wt_col = matrix(1, 10, 10)),
                   "alternating")
})


test_that("wt_row preserves row and column multisets", {

      set.seed(123)
      x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

      coords <- matrix(runif(20), ncol = 2)
      d <- as.matrix(dist(coords))
      W <- exp(-d / 0.3)

      for (method in c("curvecat", "swapcat", "tswapcat")) {
            result <- nullcat(x, method = method, n_iter = 2000,
                              wt_row = W)

            for (i in seq_len(nrow(x))) {
                  expect_equal(sort(x[i, ]), sort(result[i, ]),
                               info = paste(method, "row", i))
            }
            for (j in seq_len(ncol(x))) {
                  expect_equal(sort(x[, j]), sort(result[, j]),
                               info = paste(method, "col", j))
            }
      }
})


test_that("wt_row with block-diagonal acts like strata", {

      set.seed(123)
      x <- matrix(sample(1:3, 64, replace = TRUE), nrow = 8)

      W <- matrix(0, 8, 8)
      W[1:4, 1:4] <- 1
      W[5:8, 5:8] <- 1
      diag(W) <- 0

      result <- nullcat(x, method = "curvecat", n_iter = 5000,
                        wt_row = W, output = "index")

      idx_top <- as.vector(matrix(1:64, 8)[1:4, ])
      result_top <- as.vector(result[1:4, ])
      expect_true(all(result_top %in% idx_top),
                  info = "tokens from top block stayed in top block")

      idx_bot <- as.vector(matrix(1:64, 8)[5:8, ])
      result_bot <- as.vector(result[5:8, ])
      expect_true(all(result_bot %in% idx_bot),
                  info = "tokens from bottom block stayed in bottom block")
})


test_that("wt_col with horizontal swaps actually constrains column exchanges", {

      set.seed(123)
      # 8x8 matrix, two groups of 4 columns each
      x <- matrix(sample(1:3, 64, replace = TRUE), nrow = 8)

      # Block-diagonal weight matrix over columns: cols 1-4 swap with each other,
      # cols 5-8 swap with each other, no cross-group swaps
      W <- matrix(0, 8, 8)
      W[1:4, 1:4] <- 1
      W[5:8, 5:8] <- 1
      diag(W) <- 0

      result <- nullcat(x, method = "curvecat", n_iter = 5000,
                        wt_col = W, output = "index")

      # Check that tokens from cols 1-4 stayed in cols 1-4
      idx <- matrix(1:64, 8)
      idx_left <- as.vector(idx[, 1:4])
      result_left <- as.vector(result[, 1:4])
      expect_true(all(result_left %in% idx_left),
                  info = "tokens from left block stayed in left block")

      # Check that tokens from cols 5-8 stayed in cols 5-8
      idx_right <- as.vector(idx[, 5:8])
      result_right <- as.vector(result[, 5:8])
      expect_true(all(result_right %in% idx_right),
                  info = "tokens from right block stayed in right block")
})


test_that("wt_row auto-selects vertical swaps", {

      set.seed(123)
      x <- matrix(sample(1:4, 60, replace = TRUE), nrow = 6)

      W_row <- matrix(1, 6, 6)
      expect_no_error(nullcat(x, method = "curvecat", n_iter = 100, wt_row = W_row))

      W_col <- matrix(1, 10, 10)
      expect_no_error(nullcat(x, method = "curvecat", n_iter = 100, wt_col = W_col))
})


test_that("dual-margin wt_row + wt_col works with alternating", {

      set.seed(123)
      x <- matrix(sample(1:4, 60, replace = TRUE), nrow = 6)  # 6 x 10

      W_row <- matrix(runif(36), 6, 6); W_row <- W_row + t(W_row)
      W_col <- matrix(runif(100), 10, 10); W_col <- W_col + t(W_col)

      # Both margins
      result <- nullcat(x, method = "curvecat", n_iter = 2000,
                        wt_row = W_row, wt_col = W_col)

      for (i in seq_len(nrow(x))) {
            expect_equal(sort(x[i, ]), sort(result[i, ]))
      }
      for (j in seq_len(ncol(x))) {
            expect_equal(sort(x[, j]), sort(result[, j]))
      }

      # Explicit alternating also works
      expect_no_error(
            nullcat(x, method = "curvecat", n_iter = 100,
                    swaps = "alternating", wt_row = W_row, wt_col = W_col)
      )
})

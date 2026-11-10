# wt_row / wt_col are deprecated: weight matrices are converted to groups given
# by the connected components of their nonzero entries.

test_that("wt_row/wt_col still validate their input", {
      set.seed(123)
      x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

      expect_error(nullcat(x, wt_row = "bad"), "numeric matrix")
      expect_error(nullcat(x, wt_row = matrix(1, 3, 4)), "square matrix")
      expect_error(nullcat(x, wt_row = matrix(1, 5, 5)), "matching nrow")
      expect_error(nullcat(x, wt_col = matrix(1, 5, 5)), "matching ncol")

      W <- matrix(1, 10, 10); W[1, 2] <- -1
      expect_error(nullcat(x, wt_row = W), "non-negative")
      W <- matrix(1, 10, 10); W[1, 2] <- NA
      expect_error(nullcat(x, wt_row = W), "NA")

      # weights with nothing off the diagonal would freeze every row
      expect_error(nullcat(x, wt_row = diag(10)), "no positive off-diagonal")

      # can't combine with groups
      expect_error(nullcat(x, wt_row = matrix(1, 10, 10), row_groups = rep(1:2, 5)),
                   "only one of")
})

test_that("block weights are converted to row groups, with a deprecation warning", {
      set.seed(123)
      x <- matrix(sample(1:3, 80, replace = TRUE), nrow = 8)
      g <- rep(1:2, each = 4)
      W <- outer(g, g, "==") * runif(64)   # values don't matter, only the zero pattern

      expect_warning(r1 <- nullcat(x, wt_row = W, n_iter = 500, seed = 1),
                     "deprecated.*2 groups")
      r2 <- nullcat(x, row_groups = g, n_iter = 500, seed = 1)
      expect_identical(r1, r2)

      # asymmetric weights: a pair is connected if either direction is nonzero
      W2 <- matrix(0, 8, 8); W2[1, 2:4] <- 1; W2[8, 5:7] <- 1
      expect_warning(r3 <- nullcat(x, wt_row = W2, n_iter = 500, seed = 1), "2 groups")
      expect_identical(r3, r2)

      # column weights
      cg <- rep(1:2, length.out = 10)
      xc <- matrix(sample(1:3, 80, replace = TRUE), nrow = 8, ncol = 10)
      expect_warning(r4 <- nullcat(xc, wt_col = outer(cg, cg, "==") * 1, n_iter = 500, seed = 2),
                     "col_groups")
      expect_identical(r4, nullcat(xc, col_groups = cg, n_iter = 500, seed = 2))
})

test_that("connected weights warn that they have no effect", {
      set.seed(123)
      x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
      coords <- matrix(runif(20), ncol = 2)
      W <- exp(-as.matrix(dist(coords)) / 0.3)

      expect_warning(r <- nullcat(x, wt_row = W, n_iter = 500, seed = 3), "no effect")
      expect_identical(r, nullcat(x, n_iter = 500, seed = 3))
})

test_that("deprecated weights work through quantize, batch and trace functions", {
      set.seed(1)
      q <- matrix(rexp(200), 20)
      g <- rep(1:2, each = 10)
      W <- outer(g, g, "==") * 1

      expect_warning(r <- quantize(q, wt_row = W, n_iter = 500), "deprecated")
      for (k in 1:2) expect_equal(sort(q[g == k, ]), sort(r[g == k, ]))

      # batch converts once, so a single warning
      x <- matrix(sample(1:3, 200, replace = TRUE), 20)
      n_warn <- 0
      withCallingHandlers(nullcat_batch(x, n_reps = 3, wt_row = W, n_iter = 100),
                          warning = function(w) { n_warn <<- n_warn + 1; invokeRestart("muffleWarning") })
      expect_equal(n_warn, 1)

      expect_warning(tr <- trace_cat(x, n_iter = 100, n_chains = 2, wt_row = W), "deprecated")
      expect_equal(tr$n_blocks, 2)
})

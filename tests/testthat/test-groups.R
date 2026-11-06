cat_counts <- function(m, g, k) rowsum((m == k) * 1L, g)

test_that("row groups preserve within-group category counts for all methods", {
      set.seed(1)
      x <- matrix(sample(1:4, 30 * 12, TRUE), 30)
      rg <- rep(c("a", "b", "c"), each = 10)
      cg <- rep(1:2, each = 6)

      for (m in nullcat_methods()) {
            r <- nullcat(x, method = m, n_iter = 2000, row_groups = rg, col_groups = cg)
            # every block keeps its own multiset of categories
            for (a in unique(rg)) for (b in unique(cg)) {
                  expect_equal(sort(x[rg == a, cg == b]), sort(r[rg == a, cg == b]), info = m)
            }
            if (m %in% c("curvecat", "swapcat", "tswapcat")) {
                  # fixed margins, plus fixed counts per row group in every column
                  expect_equal(apply(x, 1, sort), apply(r, 1, sort), info = m)
                  expect_equal(apply(x, 2, sort), apply(r, 2, sort), info = m)
                  for (k in 1:4) expect_equal(cat_counts(x, rg, k), cat_counts(r, rg, k), info = m)
            }
      }
      # groups actually randomize
      expect_false(identical(nullcat(x, n_iter = 2000, row_groups = rg), x))
})

test_that("index output maps into the full matrix and stays within blocks", {
      set.seed(2)
      x <- matrix(sample(1:3, 20 * 8, TRUE), 20)
      rg <- rep(1:4, each = 5); cg <- rep(1:2, 4)
      idx <- nullcat(x, n_iter = 1000, output = "index", row_groups = rg, col_groups = cg)
      expect_setequal(as.vector(idx), seq_along(x))
      expect_equal(rg[row(x)[idx]], rg[row(x)])
      expect_equal(cg[col(x)[idx]], cg[col(x)])
})

test_that("singleton groups are left unchanged, and labels can be any type", {
      set.seed(3)
      x <- matrix(sample(1:3, 20 * 8, TRUE), 20)
      r <- nullcat(x, n_iter = 500, row_groups = c(1:5, rep(6, 15)))
      expect_identical(r[1:5, ], x[1:5, ])
      expect_identical(nullcat(x, n_iter = 200, row_groups = seq_len(20), col_groups = seq_len(8)), x)
      expect_identical(nullcat(x, row_groups = factor(rep(c("p", "q"), 10)), seed = 1),
                       nullcat(x, row_groups = rep(1:2, 10), seed = 1))
      # quantize keeps singleton rows' values in place, in every fixed mode that can't move them
      q <- matrix(runif(20 * 8), 20)
      g <- c(1:3, rep(4, 17))
      r <- quantize(q, fixed = "cell", n_iter = 500, row_groups = g)
      expect_identical(r[1:3, ], q[1:3, ])
})

test_that("group arguments are validated", {
      x <- matrix(sample(1:3, 40, TRUE), 10)
      expect_error(nullcat(x, row_groups = 1:3), "length")
      expect_error(nullcat(x, row_groups = c(NA, rep(1, 9))), "NA")
      expect_error(nullcat(x, row_groups = matrix(1, 10, 2)), "batch")
      expect_error(quantize(matrix(runif(40), 10), fixed = "stratum", row_groups = matrix(1, 10, 2)), "quantize_batch")
      expect_error(nullcat_batch(x, n_reps = 3, row_groups = matrix(1, 10, 2)), "one column per replicate")
      expect_error(nullcat_batch(x, n_reps = 2, row_groups = matrix(1, 9, 2)), "one row")
})

test_that("quantize with groups keeps values within blocks in every fixed mode", {
      set.seed(4)
      q <- matrix(runif(30 * 12), 30)
      rg <- rep(1:3, each = 10)
      s0 <- stratify(q, n_strata = 5)
      for (f in c("cell", "stratum", "row", "col")) {
            r <- quantize(q, fixed = f, n_iter = 2000, row_groups = rg)
            for (a in 1:3) expect_equal(sort(q[rg == a, ]), sort(r[rg == a, ]), info = f)
            # strata are defined globally, and each group keeps its stratum counts per species
            s1 <- stratify(r, breaks = attr(s0, "breaks"))
            for (k in 1:5) expect_equal(cat_counts(s0, rg, k), cat_counts(s1, rg, k), info = f)
      }
      # prep reuse
      prep <- quantize_prep(q, row_groups = rg, n_iter = 500)
      r <- quantize(prep = prep)
      for (a in 1:3) expect_equal(sort(q[rg == a, ]), sort(r[rg == a, ]))
})

test_that("batch functions use one column of a groups matrix per replicate", {
      set.seed(5)
      x <- matrix(sample(1:4, 30 * 6, TRUE), 30)
      q <- matrix(runif(30 * 6), 30)
      G <- matrix(sample(1:3, 30 * 4, TRUE), 30, 4)

      b <- nullcat_batch(x, n_reps = 4, row_groups = G, n_iter = 500)
      bq <- quantize_batch(q, n_reps = 4, row_groups = G, n_iter = 500)
      for (r in 1:4) {
            for (k in 1:4) expect_equal(cat_counts(x, G[, r], k), cat_counts(b[, , r], G[, r], k))
            for (a in unique(G[, r])) expect_equal(sort(q[G[, r] == a, ]), sort(bq[G[, r] == a, , r]))
      }
      # a vector applies to every replicate
      b2 <- nullcat_batch(x, n_reps = 2, row_groups = G[, 1], n_iter = 500)
      for (r in 1:2) for (k in 1:4) expect_equal(cat_counts(x, G[, 1], k), cat_counts(b2[, , r], G[, 1], k))
      # stat functions still work
      expect_length(nullcat_batch(x, n_reps = 4, row_groups = G, stat = sum), 4)
})

test_that("trace_cat traces the largest block", {
      set.seed(6)
      x <- matrix(sample(1:4, 30 * 6, TRUE), 30)
      tr <- trace_cat(x, n_iter = 200, n_chains = 2, row_groups = c(rep(1, 20), rep(2, 10)))
      expect_equal(unname(tr$block), c(20, 6))
      tq <- trace_cat(matrix(runif(180), 30), fun = "quantize", n_iter = 200, n_chains = 2,
                      row_groups = c(rep(1, 12), rep(2, 18)), col_groups = rep(1:2, 3))
      expect_equal(unname(tq$block), c(18, 3))
      expect_null(trace_cat(x, n_iter = 200, n_chains = 2)$block)
})

test_that("commsim functions accept groups", {
      skip_if_not_installed("vegan")
      set.seed(7)
      x <- matrix(sample(1:3, 20 * 5, TRUE), 20)
      g <- rep(1:2, each = 10)
      sims <- simulate(vegan::nullmodel(x, nullcat_commsim(n_iter = 500, row_groups = g)), nsim = 3)
      for (r in 1:3) for (k in 1:3) expect_equal(cat_counts(x, g, k), cat_counts(sims[, , r], g, k))
})

test_that("group_sizes averages group membership counts", {
      expect_equal(group_sizes(c("a", "a", "b")), c(2, 2, 1))
      G <- cbind(c(1, 1, 2), c(1, 2, 2))
      expect_equal(group_sizes(G), c(1.5, 2, 1.5))
      expect_equal(group_sizes(matrix(1, 1, 3)), 1)
})

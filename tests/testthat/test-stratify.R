test_that("stratify preserves matrix dimensions", {
      set.seed(123)
      m <- matrix(runif(100), nrow = 10)

      s <- stratify(m, n_strata = 5)
      expect_equal(dim(s), dim(m))
      expect_true(is.matrix(s))
      expect_true(all(s == floor(s)))  # Values are whole numbers

      # With user-supplied breaks
      s2 <- stratify(m, breaks = c(0.25, 0.5, 0.75))
      expect_equal(dim(s2), dim(m))
      expect_true(is.matrix(s2))
})


test_that("stratify preserves dimnames", {
      set.seed(123)
      m <- matrix(runif(100), nrow = 10,
                  dimnames = list(paste0("row", 1:10),
                                  paste0("col", 1:10)))

      s <- stratify(m, n_strata = 5)
      expect_equal(dimnames(s), dimnames(m))

      # With user-supplied breaks
      s2 <- stratify(m, breaks = c(0.25, 0.5, 0.75))
      expect_equal(dimnames(s2), dimnames(m))
})


test_that("stratify works with vectors", {
      x <- c(0.1, 0.5, 1.2, 3.4, 5.6, 10.2)

      s <- stratify(x, n_strata = 3)
      expect_length(s, length(x))
      expect_equal(min(s), 1L)
      expect_equal(max(s), 3L)
})


test_that("stratify creates correct number of strata", {
      set.seed(123)
      x <- runif(1000)

      for (n in 2:10) {
            s <- stratify(x, n_strata = n)
            expect_equal(length(unique(s)), n)
      }
})


test_that("stratify with user-supplied breaks", {
      x <- c(0.1, 0.3, 0.5, 0.7, 0.9)

      s <- stratify(x, breaks = c(0.5))
      expect_equal(as.integer(s), c(1L, 1L, 1L, 2L, 2L))

      s2 <- stratify(x, breaks = c(0.25, 0.5, 0.75))
      expect_equal(length(unique(s2)), 4)

      # Check breaks attribute
      breaks <- attr(s2, "breaks")
      expect_true(all(c(-Inf, 0.25, 0.5, 0.75, Inf) %in% breaks))
})


test_that("stratify with zero_stratum = TRUE", {
      x <- c(0, 0, 0.1, 0.5, 1.2, 3.4, 5.6, 10.2)

      s <- stratify(x, n_strata = 4, zero_stratum = TRUE)

      # Zeros should be in stratum 1
      expect_equal(s[x == 0], rep(1L, sum(x == 0)))

      # Non-zeros should be in strata 2-4
      expect_true(all(s[x > 0] >= 2L))
      expect_true(all(s[x > 0] <= 4L))

      # Should have all 4 strata
      expect_equal(length(unique(s)), 4)

      # Check breaks include 0
      breaks <- attr(s, "breaks")
      expect_true(0 %in% breaks)
})


test_that("stratify with zero_stratum and user breaks", {
      x <- c(0, 0, 0.5, 1.0, 1.5, 2.0)

      s <- stratify(x, breaks = c(1.0), zero_stratum = TRUE)

      # Should have 3 strata: [0], (0, 1], (1, Inf]
      expect_equal(length(unique(s)), 3)
      expect_equal(s[x == 0], rep(1L, 2))

      breaks <- attr(s, "breaks")
      expect_true(0 %in% breaks)
      expect_true(1.0 %in% breaks)
})


test_that("stratify handles all-zero case", {
      x <- rep(0, 10)

      s <- stratify(x, n_strata = 3, zero_stratum = TRUE)

      expect_true(all(s == 1L))
      expect_length(s, 10)

      breaks <- attr(s, "breaks")
      expect_equal(breaks, c(-Inf, 0, Inf))
})


test_that("stratify with transform", {
      set.seed(123)
      x <- rexp(100, rate = 0.5)

      # Log transform
      s_log <- stratify(x, n_strata = 5, transform = log)
      expect_equal(length(unique(s_log)), 5)

      # Square root transform
      s_sqrt <- stratify(x, n_strata = 5, transform = sqrt)
      expect_equal(length(unique(s_sqrt)), 5)

      # Rank transform (equal occupancy)
      s_rank <- stratify(x, n_strata = 5, transform = rank)
      expect_equal(length(unique(s_rank)), 5)

      # Rank should create more equal-sized strata
      tab_rank <- table(s_rank)
      expect_true(max(tab_rank) - min(tab_rank) <= 1)
})


test_that("stratify with transform and zero_stratum", {
      x <- c(0, 0, 0.1, 0.5, 1.2, 3.4, 5.6, 10.2)

      s <- stratify(x, n_strata = 4, transform = log1p, zero_stratum = TRUE)

      # Zeros should still be in stratum 1
      expect_equal(s[x == 0], rep(1L, 2))

      # Should have 4 strata total
      expect_equal(length(unique(s)), 4)
})


test_that("stratify with offset", {
      set.seed(123)
      x <- runif(100)

      s1 <- stratify(x, n_strata = 5, offset = 0)
      s2 <- stratify(x, n_strata = 5, offset = 0.5)
      s3 <- stratify(x, n_strata = 5, offset = -0.5)

      # Different offsets should produce different stratifications
      expect_false(identical(s1, s2))
      expect_false(identical(s1, s3))
      expect_false(identical(s2, s3))

      # But all should have 5 strata
      expect_equal(length(unique(s1)), 5)
      expect_equal(length(unique(s2)), 5)
      expect_equal(length(unique(s3)), 5)
})


test_that("stratify offset validation", {
      x <- runif(10)

      expect_error(stratify(x, offset = 1.5), "between -1 and 1")
      expect_error(stratify(x, offset = -1.5), "between -1 and 1")

      # Boundary values should work
      expect_no_error(stratify(x, offset = 1))
      expect_no_error(stratify(x, offset = -1))
})


test_that("stratify returns breaks attribute", {
      x <- runif(100)

      s <- stratify(x, n_strata = 5)
      breaks <- attr(s, "breaks")

      expect_true(!is.null(breaks))
      expect_true(is.numeric(breaks))
      expect_equal(length(breaks), 6)  # n_strata + 1
      expect_equal(breaks[1], -Inf)
      expect_equal(breaks[length(breaks)], Inf)
      expect_true(all(diff(breaks) > 0))  # monotonically increasing
})


test_that("stratify breaks are in original space when transformed", {
      set.seed(123)
      x <- rexp(100)

      s <- stratify(x, n_strata = 5, transform = log)
      breaks <- attr(s, "breaks")

      # Breaks should be in original (exp) space, not log space
      # Each break should correspond to the max value in that stratum
      for (i in 1:5) {
            stratum_values <- x[s == i]
            if (i < 5) {
                  # All values in stratum i should be <= breaks[i+1]
                  expect_true(all(stratum_values <= breaks[i + 1]))
            }
      }
})


test_that("stratify handles edge cases", {
      # Single value - can only create 1 stratum
      x <- rep(1, 10)
      s <- stratify(x, n_strata = 3)
      expect_true(all(s == 1L))
      expect_equal(length(unique(s)), 1)

      # Two distinct values - can create at most 2 strata
      x <- c(rep(1, 5), rep(2, 5))
      s <- stratify(x, n_strata = 5)
      expect_true(length(unique(s)) <= 2)

      # Very small values
      x <- runif(100, 0, 1e-10)
      expect_no_error(s <- stratify(x, n_strata = 3))
      expect_equal(length(unique(s)), 3)

      # Very large values
      x <- runif(100, 1e10, 1e11)
      expect_no_error(s <- stratify(x, n_strata = 5))
      expect_equal(length(unique(s)), 5)
})

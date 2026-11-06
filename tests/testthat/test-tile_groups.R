test_that("tile_groups returns one tiling per column", {
      set.seed(1)
      xy <- cbind(runif(100, 0, 10), runif(100, 0, 10))
      g <- tile_groups(xy, spacing = 3, n = 7, lonlat = FALSE)
      expect_true(is.integer(g))
      expect_equal(dim(g), c(100, 7))
      expect_false(anyNA(g))
      # labels recoded to 1..k within each column
      for (r in 1:7) expect_equal(sort(unique(g[, r])), seq_len(max(g[, r])))
      # one spacing per tiling
      g2 <- tile_groups(xy, spacing = c(1, 100), n = 2, lonlat = FALSE)
      expect_gt(max(g2[, 1]), max(g2[, 2]))
})

test_that("planar tiles are hexagons of the requested spacing", {
      set.seed(2)
      xy <- cbind(runif(400, 0, 50), runif(400, 0, 50))
      s <- 10
      g <- tile_groups(xy, spacing = s, n = 100, lonlat = FALSE)
      D <- as.matrix(dist(xy))
      P <- Reduce(`+`, lapply(1:100, function(r) outer(g[, r], g[, r], "=="))) / 100
      ut <- upper.tri(D)
      # nobody shares a tile beyond the hexagon diameter
      expect_lte(max(D[ut][P[ut] > 0]), 2 * s / sqrt(3) + 1e-8)
      # co-tile probability ~0.5 at ~0.42 * spacing
      near <- abs(D[ut] - 0.42 * s) < 0.5
      expect_equal(mean(P[ut][near]), 0.5, tolerance = 0.06)
})

test_that("geodesic tiles assign each site to its nearest lattice point", {
      set.seed(3)
      for (N in c(20, 500, 50000)) {
            ll <- cbind(runif(300, -180, 180), asin(runif(300, -1, 1)) * 180 / pi)
            P <- nullcat:::lonlat_to_unit(ll) %*% nullcat:::random_rotation()
            k <- nullcat:::fib_nearest(P, N)
            S <- nullcat:::fib_point(0:(N - 1), N)
            kb <- max.col(tcrossprod(P, S), ties.method = "first") - 1
            expect_equal(rowSums((S[k + 1, ] - P)^2), rowSums((S[kb + 1, ] - P)^2), tolerance = 1e-12)
      }
})

test_that("geodesic tiles match the requested spacing", {
      set.seed(4)
      ll <- cbind(runif(500, -5, 5), runif(500, -5, 5))
      s <- 300e3
      g <- tile_groups(ll, spacing = s, n = 100)
      V <- nullcat:::lonlat_to_unit(ll)
      D <- 6371008.8 * acos(pmin(tcrossprod(V), 1))
      P <- Reduce(`+`, lapply(1:100, function(r) outer(g[, r], g[, r], "=="))) / 100
      ut <- upper.tri(D)
      near <- abs(D[ut] - 0.42 * s) < 15e3
      expect_equal(mean(P[ut][near]), 0.5, tolerance = 0.06)
      expect_lte(max(D[ut][P[ut] > 0]), 1.4 * s)
})

test_that("tile_groups validates input", {
      xy <- cbind(1:5, 1:5)
      expect_error(tile_groups(xy[, 1], 1, lonlat = FALSE), "two-column")
      expect_error(tile_groups(xy, -1, lonlat = FALSE), "positive")
      expect_error(tile_groups(xy, c(1, 2, 3), n = 2, lonlat = FALSE), "spacing")
      expect_error(tile_groups(xy, 1, n = 0, lonlat = FALSE), "positive integer")
      expect_error(tile_groups(cbind(0, 100), 1e5), "lonlat = FALSE")
      expect_error(tile_groups(rbind(c(NA, 1)), 1, lonlat = FALSE), "missing")
})

test_that("tiles plug into batch randomization", {
      set.seed(5)
      xy <- cbind(runif(60, 0, 10), runif(60, 0, 10))
      x <- matrix(sample(1:3, 60 * 5, TRUE), 60)
      g <- tile_groups(xy, spacing = 4, n = 5, lonlat = FALSE)
      sims <- nullcat_batch(x, n_reps = 5, row_groups = g, n_iter = 300)
      for (r in 1:5) for (k in 1:3) {
            expect_equal(rowsum((x == k) * 1L, g[, r]), rowsum((sims[, , r] == k) * 1L, g[, r]))
      }
})

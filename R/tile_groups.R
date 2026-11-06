#' Random hexagonal tilings for spatially constrained null models
#'
#' Generates random groupings of sites by tiling space with regular hexagons
#' that are randomly positioned and oriented, independently for each column of
#' the result. Passing the result as `row_groups` to [nullcat_batch()] or
#' [quantize_batch()] randomizes each replicate within the tiles of its own
#' tiling, so that sites only exchange values with nearby sites, and the
#' resulting null distribution blends over many tile placements.
#'
#' @param coords A two-column numeric matrix (or data frame) of site
#'   coordinates, one row per row of the community matrix. With
#'   `lonlat = TRUE`, longitude and latitude in decimal degrees; otherwise
#'   planar (projected) x and y coordinates.
#' @param spacing Distance between the centres of neighbouring tiles. With
#'   `lonlat = TRUE` this is in metres (great-circle distance on a spherical
#'   Earth of radius 6,371,008.8 m); otherwise it is in the units of `coords`.
#'   Either a single value, or a vector of length `n` giving the spacing for
#'   each tiling (e.g. random draws, for a blend of scales).
#' @param n Number of tilings (columns of the result). Use one per
#'   randomization replicate.
#' @param lonlat Logical: are `coords` longitude/latitude (`TRUE`, the default)
#'   or planar coordinates (`FALSE`)?
#'
#' @details
#' Each tiling is built from a set of regularly spaced seed points, and every
#' site is assigned to its nearest seed. With `lonlat = TRUE`, seeds form a
#' Fibonacci lattice on the sphere with a uniformly random 3-D rotation for each
#' tiling. This gives near-hexagonal, near-equal-area tiles of any size,
#' anywhere on the globe, with no projection distortion; nearest seeds are
#' found by a fast exact search, so fine spacings are no more costly than
#' coarse ones. With `lonlat = FALSE`, seeds form a planar hexagonal grid with
#' a random rotation and offset for each tiling.
#'
#' Because tile placement is random, two sites a distance `d` apart share a
#' tile in some tilings and not others. Across tilings, the probability that
#' they share a tile falls from 1 at `d = 0` to about 0.5 at
#' `d = 0.42 * spacing`, and to 0 at about the tile diameter: `1.15 * spacing`
#' for a regular hexagon (the planar grid). Tiles of the spherical lattice all
#' have nearly the same area, but some are moderately elongated, so a few
#' pairs of sites up to about `1.3 * spacing` apart can occasionally share a
#' tile.
#' This "co-tile" probability acts as the distance kernel of the null model:
#' sites are pooled with their neighbours, and never with distant sites.
#' Supplying a vector of spacings blends tiles of several sizes, which smooths
#' and lengthens the kernel.
#'
#' Sites near edges, gaps, or coastlines tend to fall in tiles with few other
#' sites. A site alone in its tile cannot be changed by fixed-margin
#' randomizations in that replicate. Use [group_sizes()] to check the average
#' pool size of each site before running randomizations. When randomizing
#' within groups, `n_iter` is applied to every tile, so it only needs to be
#' large enough for the largest tile; [suggest_n_iter()] given the same groups
#' traces the largest tile.
#'
#' @return An integer matrix with one row per site and `n` columns. Each
#'   column is one tiling, giving the tile membership of every site (tile
#'   labels are arbitrary and only meaningful within a column).
#'
#' @seealso [group_sizes()], [nullcat_batch()], [quantize_batch()]
#'
#' @examples
#' # planar coordinates
#' set.seed(1)
#' xy <- cbind(runif(300, 0, 100), runif(300, 0, 100))
#' tiles <- tile_groups(xy, spacing = 30, n = 99, lonlat = FALSE)
#' summary(group_sizes(tiles))
#'
#' # randomize a community matrix within a different tiling in each replicate
#' comm <- matrix(rbinom(300 * 20, 1, 0.3), 300)
#' sims <- nullcat_batch(comm, n_reps = 99, row_groups = tiles, n_iter = 500)
#'
#' # geographic coordinates, with 200 km spacing
#' ll <- cbind(runif(300, -124, -114), runif(300, 32, 42))
#' tiles_ll <- tile_groups(ll, spacing = 200e3, n = 99)
#'
#' @export
tile_groups <- function(coords, spacing, n = 1L, lonlat = TRUE) {

      coords <- as.matrix(coords)
      if (!is.numeric(coords) || ncol(coords) != 2L) {
            stop("`coords` must be a two-column numeric matrix.", call. = FALSE)
      }
      if (nrow(coords) < 1L) stop("`coords` must have at least one row.", call. = FALSE)
      if (anyNA(coords) || any(!is.finite(coords))) {
            stop("`coords` must not contain missing or infinite values.", call. = FALSE)
      }
      n <- as.integer(n)
      if (length(n) != 1L || is.na(n) || n < 1L) stop("`n` must be a positive integer.", call. = FALSE)
      if (!is.numeric(spacing) || !(length(spacing) %in% c(1L, n)) ||
          anyNA(spacing) || any(!is.finite(spacing)) || any(spacing <= 0)) {
            stop("`spacing` must be a positive number, or a vector of `n` positive numbers.",
                 call. = FALSE)
      }
      spacing <- rep_len(spacing, n)

      if (lonlat) {
            if (any(abs(coords[, 2]) > 90) || any(abs(coords[, 1]) > 360)) {
                  stop("With `lonlat = TRUE`, `coords` must be longitude and latitude in degrees ",
                       "(use `lonlat = FALSE` for projected coordinates).", call. = FALSE)
            }
            R <- 6371008.8
            V <- lonlat_to_unit(coords)
            assign_fun <- function(s) {
                  N <- max(2, round(4 * pi * R^2 / (sqrt(3) / 2 * s^2)))
                  fib_nearest(V %*% random_rotation(), N)
            }
      } else {
            assign_fun <- function(s) hex_assign(coords, s)
      }

      out <- matrix(NA_integer_, nrow(coords), n, dimnames = list(rownames(coords), NULL))
      for (r in seq_len(n)) {
            k <- assign_fun(spacing[r])
            out[, r] <- match(k, unique(k))
      }
      out
}


# internal -- lon/lat in degrees to unit vectors
lonlat_to_unit <- function(ll) {
      lon <- ll[, 1] * pi / 180
      lat <- ll[, 2] * pi / 180
      cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}


# internal -- uniformly random rotation matrix in SO(3)
random_rotation <- function() {
      qrM <- qr(matrix(stats::rnorm(9), 3))
      Q <- qr.Q(qrM) %*% diag(sign(diag(qr.R(qrM))))
      if (det(Q) < 0) Q[, 1] <- -Q[, 1]
      Q
}


# internal -- point k (0-based) of an N-point spherical Fibonacci lattice
fib_point <- function(k, N) {
      z <- 1 - (2 * k + 1) / N
      r <- sqrt(pmax(0, 1 - z * z))
      a <- 2 * pi * ((k * ((3 - sqrt(5)) / 2)) %% 1)   # k * golden angle
      cbind(r * cos(a), r * sin(a), z)
}


# internal -- index of the nearest Fibonacci lattice point to each row of P
# (unit vectors). Starts from the lattice point at the same height and descends
# greedily using index offsets equal to Fibonacci numbers, which connect each
# lattice point to its neighbours at every scale.
fib_nearest <- function(P, N) {
      k <- pmin(pmax(round((N * (1 - P[, 3]) - 1) / 2), 0), N - 1)
      d <- rowSums((fib_point(k, N) - P)^2)
      f <- c(1, 2)
      while (f[length(f)] < N) f <- c(f, f[length(f)] + f[length(f) - 1L])
      f <- f[f < N]
      offsets <- c(rbind(f, -f))
      active <- seq_len(nrow(P))
      while (length(active)) {
            moved <- logical(length(active))
            for (o in offsets) {
                  kk <- k[active] + o
                  ok <- which(kk >= 0 & kk < N)
                  if (!length(ok)) next
                  a <- active[ok]
                  dd <- rowSums((fib_point(kk[ok], N) - P[a, , drop = FALSE])^2)
                  better <- which(dd < d[a])
                  if (length(better)) {
                        k[a[better]] <- kk[ok][better]
                        d[a[better]] <- dd[better]
                        moved[ok[better]] <- TRUE
                  }
            }
            active <- active[moved]
      }
      k
}


# internal -- assign planar points to cells of a hexagonal grid with centre
# spacing `s`, random rotation, and random offset (cube-coordinate rounding)
hex_assign <- function(xy, s) {
      th <- stats::runif(1, 0, pi / 3)
      p <- xy %*% matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
      u <- stats::runif(2)   # uniform offset over one lattice cell
      p <- sweep(p, 2, u[1] * c(s, 0) + u[2] * c(s / 2, s * sqrt(3) / 2))
      side <- s / sqrt(3)
      qf <- (sqrt(3) / 3 * p[, 1] - p[, 2] / 3) / side
      rf <- (2 / 3 * p[, 2]) / side
      yf <- -qf - rf
      rq <- round(qf); rr <- round(rf); ry <- round(yf)
      dq <- abs(rq - qf); dr <- abs(rr - rf); dy <- abs(ry - yf)
      fix_q <- dq > dy & dq > dr
      fix_r <- !fix_q & dr >= dy
      rq[fix_q] <- -ry[fix_q] - rr[fix_q]
      rr[fix_r] <- -rq[fix_r] - ry[fix_r]
      paste(rq, rr)
}

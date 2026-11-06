#' Trial-swap categorical randomization (tswapcat)
#'
#' The trial-swap ("tswap") algorithm is a fixed-fixed randomization that repeatedly
#' attempts random 2x2 swaps until a valid one is found in each iteration,
#' reducing the number of wasted draws compared to the simple swap.
#' `tswapcat()` extends this logic to categorical matrices.
#'
#' Unlike [curvecat()] and [swapcat()], whose symmetric moves converge to a
#' uniform draw from all matrices with the observed row and column category
#' multisets, redrawing until a valid swap is found biases sampling toward
#' matrices with many available swaps (Miklos & Podani 2004). The stationary
#' distribution of `tswapcat()` is therefore not uniform, and it is provided
#' mainly for methodological comparison.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#'
#' @references
#' Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns.
#' *Ecology*, 81(9), 2606-2621.
#'
#' Miklos, I. & Podani, J. (2004). Randomization of presence-absence matrices:
#' comments and new algorithms. *Ecology*, 85(1), 86-92.
#'
#' @seealso [curvecat()] for an efficient algorithm that samples matrices
#'   uniformly.
#'
#' @examples
#' set.seed(123)
#' x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
#' x_rand <- tswapcat(x, n_iter = 1000)
#'
#' @export
tswapcat <- function(x, n_iter = 1000L, output = c("category", "index"), swaps = "auto",
                     row_groups = NULL, col_groups = NULL,
                     wt_row = NULL, wt_col = NULL, seed = NULL) {
      nullcat(x, method = "tswapcat", n_iter = n_iter, output = output,
              swaps = swaps, row_groups = row_groups, col_groups = col_groups,
              wt_row = wt_row, wt_col = wt_col, seed = seed)
}

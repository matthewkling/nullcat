#' Trial-swap categorical randomization (tswapcat)
#'
#' The trial-swap ("tswap") algorithm is a fixed–fixed randomization that repeatedly
#' attempts random 2×2 swaps until a valid one is found in each iteration,
#' reducing the number of wasted draws compared to the simple swap.
#' \code{tswapcat()} extends this logic to categorical matrices.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#'
#' @references
#' Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns.
#' *Ecology*, 81(9), 2606–2621.
#'
#' Miklós, I. & Podani, J. (2004). Randomization of presence–absence matrices:
#' comments and new algorithms. *Ecology*, 85(1), 86–92.
#'
#' Gotelli, N. J. & Entsminger, G. L. (2003). *EcoSim: Null models software for
#' ecology* (Version 7.0). Acquired Intelligence Inc. & Kesey-Bear, Jericho (VT).
#'
#' @examples
#' set.seed(123)
#' x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
#'
#' # Randomize using swap algorithm
#' x_rand <- tswapcat(x, n_iter = 1000)
#'
#' # Verify fixed-fixed constraint (row and column margins preserved)
#' all.equal(sort(x[1, ]), sort(x_rand[1, ]))
#' all.equal(sort(x[, 1]), sort(x_rand[, 1]))
#'
#' @export
tswapcat <- function(x, n_iter = 1000L, output = c("category", "index"), seed = NULL) {
      nullcat(x, method = "tswapcat", n_iter = n_iter, output = output, seed = seed)
}

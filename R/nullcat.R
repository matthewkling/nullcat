
#' Categorical matrix randomization
#'
#' Categorical generalizations of binary community null model algorithms.
#'
#' @param x A matrix of categorical data, encoded as integers.
#'   Values should represent category or stratum membership for each cell.
#' @param method Character specifying the randomization algorithm to use.
#'   Options include the following; see details and linked functions
#'   for more info.
#'   \itemize{
#'     \item \code{"curvecat"}: categorical analog to `curveball`;
#'          see \link{curvecat} for details.
#'     \item \code{"swapcat"}: categorical analog to `swap`;
#'          see \link{swapcat} for details.
#'     \item \code{"tswapcat"}: categorical analog to `tswap`;
#'          see \link{tswapcat} for details.
#'     \item \code{"r0cat"}: categorical analog to `r0`;
#'          see \link{r0cat} for details.
#'     \item \code{"c0cat"}: categorical analog to `c0`;
#'          see \link{c0cat} for details.
#'   }
#' @param n_iter Number of iterations. Default is 1000. Larger values yield
#'   more thorough mixing. Ignored for non-sequential methods. Minimum
#'   burn-in times can be estimated with \link{suggest_n_iter}.
#' @param output Character indicating type of result to return:
#'   \itemize{
#'     \item \code{"category"} (default) returns randomized matrix
#'     \item \code{"index"} returns an index matrix describing where original
#'          entries (a.k.a. "tokens") moved. Useful mainly for testing, and for
#'          applications like \code{quantize()} that care about token tracking
#'          in addition to generic integer categories.
#'   }
#' @param swaps Character string controlling the direction of token movement.
#'   Only used when method is `curvecat`, `swapcat`, or `tswapcat`.
#'   Affects the result only when \code{output = "index"}, otherwise it only affects
#'   computation speed. Options include:
#'   \itemize{
#'     \item \code{"vertical"}: Tokens move between rows (stay within columns).
#'     \item \code{"horizontal"}: Tokens move between columns (stay within rows).
#'     \item \code{"alternating"}: Tokens move in both dimensions, alternating between
#'       vertical and horizontal swaps. Provides full 2D mixing without preserving
#'       either row or column token sets.
#'     \item \code{"auto"} (default): For \code{output = "category"},
#'       automatically selects the fastest option based on matrix dimensions. For
#'       \code{output = "index"}, defaults to \code{"alternating"} for full mixing.
#'   }
#' @param seed Integer used to seed random number generator, for reproducibility.
#'
#' @details
#' `curvecat`, `swapcat`, and `tswapcat` are sequential algorithms that hold
#' category multisets fixed in every row and column. These three algorithms
#' typically reach the same stationary distribution. They differ primarily in
#' efficiency, with `curvecat` being the most efficient (i.e. fewest steps to
#' become fully mixed); `swapcat` and `tswapcat` are thus useful mainly for
#' methodological comparison.
#'
#' The `r0cat` algorithm holds category multisets fixed in rows but not columns,
#' while `c0cat` does the opposite.
#'
#' Note that categorical null models are for cell-level categorical data. Site-level
#' attributes (e.g., land cover) or species-level attributes (e.g., functional
#' traits) should be analyzed using different approaches.
#'
#' @return A matrix of the same dimensions as \code{x}, either randomized
#'   categorical values (when \code{output = "category"}) or an integer index
#'   matrix describing the permutation of entries (when \code{output = "index"}).
#'
#' @examples
#' # Create a categorical matrix
#' set.seed(123)
#' x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
#'
#' # Randomize using curvecat method (preserves row & column margins)
#' x_rand <- nullcat(x, method = "curvecat", n_iter = 1000)
#'
#' # Check that row multisets are preserved
#' all.equal(sort(x[1, ]), sort(x_rand[1, ]))
#'
#' # Get index output showing where each cell moved
#' idx <- nullcat(x, method = "curvecat", n_iter = 1000, output = "index")
#'
#' # Use different methods
#' x_swap <- nullcat(x, method = "swapcat", n_iter = 1000)
#' x_r0 <- nullcat(x, method = "r0cat")
#'
#' # Use with a seed for reproducibility
#' x_rand1 <- nullcat(x, method = "curvecat", n_iter = 1000, seed = 42)
#' x_rand2 <- nullcat(x, method = "curvecat", n_iter = 1000, seed = 42)
#' identical(x_rand1, x_rand2)
#'
#' @export
nullcat <- function(x,
                    method = nullcat_methods(),
                    n_iter = 1000L,
                    output = c("category", "index"),
                    swaps = c("auto", "vertical", "horizontal", "alternating"),
                    seed = NULL) {

      method <- match.arg(method, NULLCAT_METHODS)
      output <- match.arg(output)
      swaps <- match.arg(swaps)
      n_iter <- as.integer(n_iter)

      x <- as.matrix(x)
      storage.mode(x) <- "integer"

      # Handle "auto" swaps setting
      if (swaps == "auto" & method %in% c("curvecat", "swapcat", "tswapcat")){
            if (output == "category") {
                  # Pick fastest based on dimensions
                  # (confirmed with benchmarking. it's likely about minimizing the inner loop
                  # iterations, not about the sampling space or cache)
                  swaps <- if (nrow(x) > ncol(x)) "vertical" else "horizontal"
            } else {
                  # For token tracking, use alternating for full 2D mixing
                  swaps <- "alternating"
            }
      }

      with_seed(seed, {
            switch(
                  method,
                  "curvecat" = curvecat_cpp(x, n_iter = n_iter, output = output, swaps = swaps),
                  "swapcat"  = swapcat_cpp(x, n_iter = n_iter, output = output, swaps = swaps),
                  "tswapcat"  = tswapcat_cpp(x, n_iter = n_iter, output = output, swaps = swaps),
                  "r0cat"  = r0cat_cpp(x, n_iter = 1, output = output),
                  "c0cat"  = c0cat_cpp(x, n_iter = 1, output = output)
            )
      })
}


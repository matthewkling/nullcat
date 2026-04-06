#' @keywords internal
NULLCAT_METHODS <- c("curvecat", "swapcat", "tswapcat", "r0cat", "c0cat")

#' Available categorical null model methods
#'
#' @return Character vector of available method names.
#' @export
nullcat_methods <- function() NULLCAT_METHODS


#' Categorical matrix randomization
#'
#' Randomize binary or categorical community matrices using categorical
#' generalizations of binary community null model algorithms. Optionally
#' constrain mixing using spatial (row) and taxonomic (column) weights.
#'
#' @param x A matrix of categorical data, encoded as integers.
#'   Values should represent category or stratum membership for each cell.
#' @param method Character specifying the randomization algorithm to use.
#'   Options include the following; see details and linked functions
#'   for more info.
#'   \itemize{
#'     \item `"curvecat"`: categorical analog to `curveball`;
#'          see [curvecat()] for details.
#'     \item `"swapcat"`: categorical analog to `swap`;
#'          see [swapcat()] for details.
#'     \item `"tswapcat"`: categorical analog to `tswap`;
#'          see [tswapcat()] for details.
#'     \item `"r0cat"`: categorical analog to `r0`;
#'          see [r0cat()] for details.
#'     \item `"c0cat"`: categorical analog to `c0`;
#'          see [c0cat()] for details.
#'   }
#' @param n_iter Number of iterations. Default is 1000. Larger values yield
#'   more thorough mixing. Ignored for non-sequential methods. Minimum
#'   burn-in times can be estimated with `suggest_n_iter()`.
#' @param output Character indicating type of result to return:
#'   \itemize{
#'     \item `"category"` (default) returns randomized matrix
#'     \item `"index"` returns an index matrix describing where original
#'          entries (a.k.a. "tokens") moved. Useful mainly for testing, and for
#'          applications like `quantize()` that care about token tracking
#'          in addition to generic integer categories.
#'   }
#' @param swaps Character string controlling the direction of token movement.
#'   Only used when method is `"curvecat"`, `"swapcat"`, or `"tswapcat"`.
#'   Affects the result only when `output = "index"`, otherwise it only affects
#'   computation speed. Options include:
#'   \itemize{
#'     \item `"vertical"`: Tokens move between rows (stay within columns).
#'     \item "`horizontal"`: Tokens move between columns (stay within rows).
#'     \item `"alternating"`: Tokens move in both dimensions, alternating between
#'       vertical and horizontal swaps. Provides full 2D mixing without preserving
#'       either row or column token sets.
#'     \item `"auto"` (default): For `output = "category"`,
#'       automatically selects the fastest option based on matrix dimensions. For
#'       `output = "index"`, defaults to `"alternating"` for full mixing.
#'       When `wt_row` or `wt_col` is supplied, defaults to the appropriate
#'       direction, or `"alternating"` if both are supplied.
#'   }
#' @param wt_row An optional square numeric matrix of non-negative weights
#'   controlling which pairs of rows are likely to exchange tokens during
#'   randomization. Must be `nrow(x)` by `nrow(x)`. This enables spatially
#'   or trait-constrained null models where nearby or similar sites exchange
#'   tokens more frequently.
#'
#'   Values are treated as relative weights (not probabilities) and are
#'   normalized internally. The diagonal is ignored. The matrix should be
#'   symmetric. Only supported for sequential methods (`curvecat`, `swapcat`,
#'   `tswapcat`).
#'
#'   When both `wt_row` and `wt_col` are supplied, `swaps` is forced to
#'   `"alternating"`, producing a Gibbs-like sweep that applies each weight
#'   matrix on its respective margin in alternation.
#' @param wt_col An optional square numeric matrix of non-negative weights
#'   controlling which pairs of columns are likely to exchange tokens during
#'   randomization. Must be `ncol(x)` by `ncol(x)`. See `wt_row` for details
#'   on weight interpretation.
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
#' traits) should be analyzed using different approaches. See vignette for details.
#'
#' @seealso [nullcat_batch()] for efficient generation of multiple randomized
#'   matrices; [nullcat_commsim()] for integration with `vegan`.
#'
#' @return A matrix of the same dimensions as `x`, either randomized
#'   categorical values (when `output = "category"`) or an integer index
#'   matrix describing the permutation of entries (when `output = "index"`).
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
#' # Spatially constrained randomization using row weights
#' coords <- cbind(runif(10), runif(10))
#' d <- as.matrix(dist(coords))
#' W <- exp(-d / 0.3)  # Gaussian distance decay
#' x_spatial <- nullcat(x, method = "curvecat", n_iter = 1000, wt_row = W)
#'
#' # Dual-margin weighting (Gibbs-like alternating)
#' W_row <- exp(-as.matrix(dist(cbind(runif(10), runif(10)))) / 0.3)
#' W_col <- exp(-as.matrix(dist(cbind(runif(10), runif(10)))) / 0.3)
#' x_dual <- nullcat(x, method = "curvecat", n_iter = 1000,
#'                   wt_row = W_row, wt_col = W_col)
#'
#' @export
nullcat <- function(x,
                    method = nullcat_methods(),
                    n_iter = 1000L,
                    output = c("category", "index"),
                    swaps = c("auto", "vertical", "horizontal", "alternating"),
                    wt_row = NULL,
                    wt_col = NULL,
                    seed = NULL) {

      method <- match.arg(method, NULLCAT_METHODS)
      output <- match.arg(output)
      swaps <- match.arg(swaps)
      n_iter <- as.integer(n_iter)

      x <- as.matrix(x)
      storage.mode(x) <- "integer"

      has_wt <- !is.null(wt_row) || !is.null(wt_col)

      if (has_wt) {
            if (!method %in% c("curvecat", "swapcat", "tswapcat")) {
                  stop("wt_row/wt_col are only supported for sequential methods ",
                       "(curvecat, swapcat, tswapcat).")
            }
      }

      # Validate wt_row
      if (!is.null(wt_row)) {
            if (!is.matrix(wt_row) || !is.numeric(wt_row)) stop("wt_row must be a numeric matrix.")
            if (nrow(wt_row) != ncol(wt_row)) stop("wt_row must be a square matrix.")
            if (nrow(wt_row) != nrow(x)) stop("wt_row must be ", nrow(x), " x ", nrow(x), " (matching nrow(x)).")
            if (any(wt_row < 0, na.rm = TRUE)) stop("wt_row must contain only non-negative values.")
            if (anyNA(wt_row)) stop("wt_row must not contain NA values.")
      }

      # Validate wt_col
      if (!is.null(wt_col)) {
            if (!is.matrix(wt_col) || !is.numeric(wt_col)) stop("wt_col must be a numeric matrix.")
            if (nrow(wt_col) != ncol(wt_col)) stop("wt_col must be a square matrix.")
            if (nrow(wt_col) != ncol(x)) stop("wt_col must be ", ncol(x), " x ", ncol(x), " (matching ncol(x)).")
            if (any(wt_col < 0, na.rm = TRUE)) stop("wt_col must contain only non-negative values.")
            if (anyNA(wt_col)) stop("wt_col must not contain NA values.")
      }

      # Resolve swaps direction based on weights
      if (swaps == "auto" & has_wt) {
            if (!is.null(wt_row) && !is.null(wt_col)) {
                  swaps <- "alternating"
            } else if (!is.null(wt_row)) {
                  swaps <- "vertical"
            } else {
                  swaps <- "horizontal"
            }
      }

      # Check consistency of weights with swaps direction
      if (!is.null(wt_row) && swaps == "horizontal") {
            stop("wt_row cannot be used with swaps = 'horizontal' ",
                 "(row weights require vertical or alternating swaps).")
      }
      if (!is.null(wt_col) && swaps == "vertical") {
            stop("wt_col cannot be used with swaps = 'vertical' ",
                 "(column weights require horizontal or alternating swaps).")
      }
      if (!is.null(wt_row) && !is.null(wt_col) && swaps != "alternating") {
            stop("Both wt_row and wt_col supplied; swaps must be 'alternating'.")
      }

      # Handle "auto" swaps setting (when no weights supplied)
      if (swaps == "auto" & method %in% c("curvecat", "swapcat", "tswapcat")){
            if (output == "category") {
                  swaps <- if (nrow(x) > ncol(x)) "vertical" else "horizontal"
            } else {
                  swaps <- "alternating"
            }
      }

      with_seed(seed, {
            switch(
                  method,
                  "curvecat" = curvecat_cpp(x, n_iter = n_iter, output = output,
                                            swaps = swaps, wt_row = wt_row,
                                            wt_col = wt_col),
                  "swapcat"  = swapcat_cpp(x, n_iter = n_iter, output = output,
                                           swaps = swaps, wt_row = wt_row,
                                           wt_col = wt_col),
                  "tswapcat"  = tswapcat_cpp(x, n_iter = n_iter, output = output,
                                             swaps = swaps, wt_row = wt_row,
                                             wt_col = wt_col),
                  "r0cat"  = r0cat_cpp(x, n_iter = 1, output = output),
                  "c0cat"  = c0cat_cpp(x, n_iter = 1, output = output)
            )
      })
}

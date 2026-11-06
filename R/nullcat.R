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
#' constrain randomization to occur within groups of rows (e.g. regions) and/or
#' columns (e.g. clades).
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
#'   burn-in times can be estimated with `suggest_n_iter()`. When `row_groups`
#'   or `col_groups` is used, `n_iter` iterations are run separately on every
#'   block, so it should be sized for a block rather than for the full matrix.
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
#'   }
#' @param row_groups,col_groups Optional vectors assigning each row (for
#'   `row_groups`) or column (for `col_groups`) of `x` to a group. Labels can be
#'   integers, characters, or factor levels. When supplied, the matrix is split
#'   into blocks (every combination of a row group and a column group) and each
#'   block is randomized independently with `method`, so values are only ever
#'   exchanged within a block. For algorithms that fix row and column margins,
#'   this additionally fixes the category counts of every row group within each
#'   column (and every column group within each row): for example, grouping
#'   sites into regions gives a null model in which each region keeps its own
#'   species pool. Rows or columns that are alone in their group cannot be
#'   changed by fixed-margin algorithms. To randomize within a different
#'   grouping in each replicate (e.g. random spatial tilings from
#'   [tile_groups()]), pass a groups matrix to [nullcat_batch()] or
#'   [quantize_batch()].
#' @param wt_row,wt_col Deprecated. Experimental pairwise weight matrices previously used to
#'   bias which pairs of rows or columns exchange values. But weighting pair
#'   selection changes only how quickly a randomization mixes, not the null
#'   distribution it converges to; the only lasting effect of a weight matrix
#'   is to prevent exchanges between rows (or columns) that are not connected,
#'   directly or through intermediaries, by nonzero weights. Weight matrices
#'   are therefore now converted to `row_groups` / `col_groups` given by the
#'   connected components of their nonzero entries, with a warning. Use
#'   `row_groups` / `col_groups` instead.
#' @param seed Integer used to seed random number generator, for reproducibility.
#'
#' @details
#' `curvecat`, `swapcat`, and `tswapcat` are sequential algorithms that hold
#' category multisets fixed in every row and column. `curvecat` and `swapcat`
#' propose symmetric moves, so with enough iterations they converge to the same
#' null distribution: a uniform draw from all matrices with the observed row and
#' column category multisets. `curvecat` is the most efficient of the two
#' (i.e. fewest steps to become fully mixed). `tswapcat` repeatedly redraws
#' until it finds a valid swap; this makes each iteration productive, but
#' biases the sampling toward matrices with many available swaps, so its
#' stationary distribution is not uniform (cf. Miklós & Podani 2004).
#' `swapcat` and `tswapcat` are thus useful mainly for methodological
#' comparison.
#'
#' The `r0cat` algorithm holds category multisets fixed in rows but not columns,
#' while `c0cat` does the opposite.
#'
#' Note that categorical null models are for cell-level categorical data. Site-level
#' attributes (e.g., land cover) or species-level attributes (e.g., functional
#' traits) should be analyzed using different approaches. See vignette for details.
#'
#' @references
#' Miklós, I., & Podani, J. (2004). Randomization of presence–absence matrices:
#' comments and new algorithms. *Ecology*, 85(1), 86–92.
#'
#' @seealso [nullcat_batch()] for efficient generation of multiple randomized
#'   matrices; [nullcat_commsim()] for integration with `vegan`;
#'   [tile_groups()] for random spatial groupings.
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
#' # Randomize within two regions of five sites each
#' region <- rep(c("north", "south"), each = 5)
#' x_reg <- nullcat(x, method = "curvecat", n_iter = 1000, row_groups = region)
#'
#' # Each region keeps its own category counts in every column
#' all.equal(rowsum((x == 1) * 1, region), rowsum((x_reg == 1) * 1, region))
#'
#' @export
nullcat <- function(x,
                    method = nullcat_methods(),
                    n_iter = 1000L,
                    output = c("category", "index"),
                    swaps = c("auto", "vertical", "horizontal", "alternating"),
                    row_groups = NULL,
                    col_groups = NULL,
                    wt_row = NULL,
                    wt_col = NULL,
                    seed = NULL) {

      method <- match.arg(method, NULLCAT_METHODS)
      output <- match.arg(output)
      swaps <- match.arg(swaps)
      n_iter <- as.integer(n_iter)

      x <- as.matrix(x)
      storage.mode(x) <- "integer"

      row_groups <- resolve_groups(row_groups, wt_row, nrow(x), "row_groups", "wt_row")
      col_groups <- resolve_groups(col_groups, wt_col, ncol(x), "col_groups", "wt_col")

      run <- function(m) nullcat_engine(m, method = method, n_iter = n_iter,
                                        output = output, swaps = swaps)

      with_seed(seed, {
            if (is.null(row_groups) && is.null(col_groups)) {
                  run(x)
            } else {
                  apply_blocks(x, row_groups, col_groups, run, output = output)
            }
      })
}


# internal -- run one randomization algorithm on a whole matrix.
nullcat_engine <- function(x, method, n_iter, output, swaps) {

      # Sequential algorithms can't change a matrix with a single row or column
      # (or with no iterations). Handle that here: the C++ engines return the
      # input categories in these cases, even when an index is requested.
      if (method %in% c("curvecat", "swapcat", "tswapcat") &&
          (nrow(x) < 2L || ncol(x) < 2L || n_iter < 1L)) {
            if (output == "index") return(matrix(seq_along(x), nrow(x), ncol(x)))
            return(x)
      }

      # Handle "auto" swaps setting
      if (swaps == "auto" & method %in% c("curvecat", "swapcat", "tswapcat")){
            if (output == "category") {
                  swaps <- if (nrow(x) > ncol(x)) "vertical" else "horizontal"
            } else {
                  swaps <- "alternating"
            }
      }

      switch(
            method,
            "curvecat" = curvecat_cpp(x, n_iter = n_iter, output = output, swaps = swaps),
            "swapcat"  = swapcat_cpp(x, n_iter = n_iter, output = output, swaps = swaps),
            "tswapcat" = tswapcat_cpp(x, n_iter = n_iter, output = output, swaps = swaps),
            "r0cat"    = r0cat_cpp(x, n_iter = 1, output = output),
            "c0cat"    = c0cat_cpp(x, n_iter = 1, output = output)
      )
}

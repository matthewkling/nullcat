#' Generate a batch of null matrices using quantize()
#'
#' Runs the stratified null model implemented in [quantize()] repeatedly,
#' generating a batch of randomized matrices or, optionally, a batch of
#' summary statistics computed from those matrices.
#'
#' @param x Community matrix (species × sites, or any numeric matrix).
#' @param n_reps Number of randomizations to generate. Default is `999`.
#' @param stat Optional summary function taking a matrix and returning a numeric
#'        statistic (e.g. `rowSums` with abundance data would give total abundance per site).
#'        If `NULL` (default), the function returns the full set of randomized matrices.
#' @param n_cores Number of compute cores to use for parallel processing. Default is `1`.
#' @param ... Additional arguments passed to `quantize()`,
#'        (e.g. `method`, `breaks`, `n_strata`, `transform`, `offset`, `zero_stratum`,
#'        `fixed`, `n_iter`, etc.). `row_groups` and `col_groups` may be given
#'        either as vectors (the same grouping for every replicate) or as
#'        matrices with one column per replicate (e.g. random tilings from
#'        [tile_groups()]), in which case replicate `r` uses column `r`.
#' @param seed Integer used to seed random number generator, for reproducibility.
#'
#' @return If `stat` is `NULL`, returns a 3D array (rows × cols × n_reps).
#'   If `stat` is not `NULL`, returns a numeric array of statistic values
#'   (dimensionality depends on `stat`).
#'
#' @examples
#' set.seed(123)
#' x <- matrix(runif(100), nrow = 10)
#'
#' # Generate 99 randomized matrices
#' nulls <- quantize_batch(x, n_reps = 99, method = "curvecat", n_iter = 100)
#'
#' # Or compute a statistic on each
#' row_sums <- nullcat_batch(x, n_reps = 99, stat = rowSums,
#'                           method = "curvecat", n_iter = 100)
#'
#' # Specify multiple cores for parallel processing
#' nulls <- quantize_batch(x, n_reps = 99, n_iter = 100, n_cores = 2)
#'
#' @export
quantize_batch <- function(x, n_reps = 999L, stat = NULL, n_cores = 1L, seed = NULL, ...) {

      # one-time overhead
      prep <- quantize_prep(as.matrix(x), ...)

      check_groups_reps(prep$row_groups, n_reps, "row_groups")
      check_groups_reps(prep$col_groups, n_reps, "col_groups")

      # define per-rep function
      one <- function(i) quantize_draw(prep, groups_rep(prep$row_groups, i),
                                       groups_rep(prep$col_groups, i))
      if (is.null(stat)) {
            fun <- one
      } else {
            fun <- function(i) stat(one(i))
      }

      with_seed(seed, {
            sims <- mc_replicate(n_reps, fun, n_cores = n_cores)
      })

      sims
}

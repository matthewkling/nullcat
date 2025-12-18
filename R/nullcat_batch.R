#' Generate a batch of null matrices using nullcat()
#'
#' Runs the categorical null model implemented in [nullcat()] repeatedly,
#' generating a batch of randomized matrices or, optionally, a batch of
#' summary statistics computed from those matrices. This is the categorical
#' analog of `quantize_batch()`.
#'
#' @param x Community matrix (sites × species) or any categorical matrix of integers.
#' @param n_reps Number of randomizations to generate. Default is `999`.
#' @param stat Optional summary function taking a matrix and returning a numeric
#'        statistic. If `NULL` (default), the function returns the full set of
#'        randomized matrices.
#' @param n_cores Number of compute cores to use for parallel processing. Default is `1`.
#' @param seed Integer used to seed random number generator, for reproducibility.
#' @param ... Additional arguments passed to [nullcat()]
#'        (e.g. `method`, `n_iter`, `output`).
#'
#' @return If `stat` is `NULL`, returns a 3D array (rows × cols × n_reps).
#'   If `stat` is not `NULL`, returns a numeric array of statistic values
#'   (dimensionality depends on `stat`).
#'
#' @examples
#' set.seed(123)
#' x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
#'
#' # Generate 99 randomized matrices
#' nulls <- nullcat_batch(x, n_reps = 99, method = "curvecat", n_iter = 100)
#'
#' # Or compute a statistic on each
#' row_sums <- nullcat_batch(x, n_reps = 99, stat = rowSums,
#'                           method = "curvecat", n_iter = 100)
#'
#' # Specify multiple cores for parallel processing
#' nulls <- nullcat_batch(x, n_reps = 99, n_iter = 100, n_cores = 2)
#'
#' @export
nullcat_batch <- function(x, n_reps = 999L, stat = NULL, n_cores = 1L, seed = NULL, ...) {

      if (is.null(stat)) {
            fun <- function() nullcat(x, ...)
      } else {
            fun <- function() stat(nullcat(x, ...))
      }

      with_seed(seed, {
            sims <- mc_replicate(n_reps, fun, n_cores = n_cores)
      })

      sims
}

#' Categorical curveball randomization (curvecat)
#'
#' Categorical generalization of the binary curveball algorithm
#' (Strona et al. 2014) to matrices of categorical data. This function
#' is a convenience wrapper around [nullcat()] with
#' `method = "curvecat"`.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#'
#' @details
#' The curvecat algorithm randomizes a categorical matrix while keeping
#' the category multisets of each row and column fixed. In other words,
#' the permuted matrix has the same set of integer values in every row
#' and every column as the original matrix, but they are permuted. It
#' operates on pairs of rows at a time, grouping differing entries by
#' unordered category pairs and redistributing the orientation of those
#' pairs while preservingn the multiset of categories within each row.
#' When there are only two categories, `curvecat()` reduces to the
#' behavior of the original binary curveball algorithm (Strona et al. 2014)
#' applied to a 0/1 matrix.
#'
#' @seealso [nullcat()], [nullcat_methods()]
#'
#' @references
#' Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz, J.
#' (2014). A fast and unbiased procedure to randomize ecological binary
#' matrices with fixed row and column totals. \emph{Nature Communications}, 5,
#' 4114.
#'
#' @examples
#' # Create a categorical matrix
#' set.seed(123)
#' x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
#'
#' # Randomize preserving row and column category multisets
#' x_rand <- curvecat(x, n_iter = 1000)
#'
#' # Verify margins are preserved
#' all.equal(sort(x[1, ]), sort(x_rand[1, ])) # row multisets preserved
#' all.equal(sort(x[, 1]), sort(x_rand[, 1])) # column multisets preserved
#'
#' # Use with a seed for reproducibility
#' x_rand1 <- curvecat(x, n_iter = 1000, seed = 42)
#' x_rand2 <- curvecat(x, n_iter = 1000, seed = 42)
#' identical(x_rand1, x_rand2)
#'
#' @export
curvecat <- function(x, n_iter = 1000L, output = "category", swaps = "auto", seed = NULL) {
      nullcat(x, method = "curvecat", n_iter = n_iter, output = output, swaps = swaps, seed = seed)
}

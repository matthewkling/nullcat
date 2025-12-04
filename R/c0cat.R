#' Column-constrained categorical randomization (c0cat)
#'
#' `c0cat()` preserves the multiset of categories within each column but
#' randomizes their positions across rows, leaving row margins free.
#' This is the categorical analog to vegan's `c0` algorithm. It is a
#' non-sequential method.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#' @examples
#' set.seed(123)
#' x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
#'
#' # Randomize within columns (column margins fixed, row margins free)
#' x_rand <- c0cat(x)
#'
#' # Verify columns are preserved but rows are not
#' all.equal(sort(x[, 1]), sort(x_rand[, 1]))
#' any(sort(x[1, ]) != sort(x_rand[1, ]))
#'
#' @export
c0cat <- function(x, n_iter = 1L, output = c("category", "index"), seed = NULL) {
      nullcat(x, method = "c0cat", n_iter = n_iter, output = output, seed = seed)
}

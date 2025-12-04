#' Row-constrained categorical randomization (r0cat)
#'
#' `r0cat()` preserves the multiset of categories within each row but
#' randomizes their positions across columns, leaving column margins free.
#' This is the categorical analog to vegan's `r0` algorithm. It is a
#' non-sequential method.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#' @examples
#' set.seed(123)
#' x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
#'
#' # Randomize within rows (row margins fixed, column margins free)
#' x_rand <- r0cat(x)
#'
#' # Verify rows are preserved but columns are not
#' all.equal(sort(x[1, ]), sort(x_rand[1, ]))
#' any(sort(x[, 1]) != sort(x_rand[, 1]))
#'
#' @export
r0cat <- function(x, n_iter = 1L, output = c("category", "index"), seed = NULL) {
      nullcat(x, method = "r0cat", n_iter = n_iter, output = output, seed = seed)
}

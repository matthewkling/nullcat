#' Row-constrained categorical randomization (r0cat)
#'
#' `r0cat()` preserves the multiset of categories within each row but
#' randomizes their positions across columns, leaving column margins free.
#' This is the categorical analog to vegan's `r0` algorithm.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#' @export
r0cat <- function(x, n_iter = 1L, output = c("category", "index")) {
      nullcat(x, method = "r0cat", n_iter = n_iter, output = output)
}

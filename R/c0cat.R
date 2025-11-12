#' Column-constrained categorical randomization (c0cat)
#'
#' `r0cat()` preserves the multiset of categories within each column but
#' randomizes their positions across rows, leaving row margins free.
#' This is the categorical analog to vegan's `c0` algorithm. It is a
#' non-sequential method.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#' @export
c0cat <- function(x, n_iter = 1L, output = c("category", "index")) {
      nullcat(x, method = "c0cat", n_iter = n_iter, output = output)
}

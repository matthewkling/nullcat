#' Categorical curveball randomization (curvecat)
#'
#' Categorical generalization of the binary curveball algorithm
#' (Strona et al.) to matrices of categorical data. This function
#' is a convenience wrapper around [nullcat()] with
#' `method = "curvecat"`.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#'
#' @details
#' The curvecat algorithm operates on pairs of rows at a time,
#' grouping differing entries by unordered category pairs and
#' redistributing the orientation of those pairs while preserving
#' the multiset of categories within each row. When there are only
#' two categories, `curvecat()` reduces to the behavior of the
#' original binary curveball algorithm applied to a 0/1 matrix.
#'
#'
#' @seealso [nullcat()], [nullcat_methods()]
#'
#' @references
#' Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz, J.
#' (2014). A fast and unbiased procedure to randomize ecological binary
#' matrices with fixed row and column totals. \emph{Nature Communications}, 5,
#' 4114.
#'
#' @export
curvecat <- function(x, n_iter = 1000L, output = c("category", "index"), seed = NULL) {
      nullcat(x, method = "curvecat", n_iter = n_iter, output = output, seed = seed)
}

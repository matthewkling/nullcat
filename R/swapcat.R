#' Categorical swap randomization (swapcat)
#'
#' Categorical generalization of the binary 2x2 swap algorithm to
#' matrices of categorical data. This function is a convenience wrapper
#' around [nullcat()] with `method = "swapcat"`.
#'
#' @inheritParams nullcat
#' @inherit nullcat return
#'
#' @details
#' The swapcat algorithm attempts random 2x2 swaps of the form:
#'
#' ```
#' a b        b a
#' b a   <->  a b
#' ```
#'
#' where \eqn{a} and \eqn{b} are distinct categories. These swaps
#' preserve the multiset of categories in each row and column.
#' With only two categories present, `swapcat()` reduces to the
#' behavior of the standard binary swap algorithm.
#'
#' @seealso [curvecat()] for an algorithm that produces equivalent results with
#'   better computational efficiency.
#'
#' @references
#' Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns.
#' *Ecology*, 81(9), 2606–2621.
#'
#' See also Gotelli & Entsminger (2003) *EcoSim: Null models software for ecology*
#' (Version 7.0) for implementation details of the binary swap algorithm.
#'
#' @examples
#' set.seed(123)
#' x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
#'
#' # Randomize using swap algorithm
#' x_rand <- swapcat(x, n_iter = 1000)
#'
#' # Verify fixed-fixed constraint (row and column margins preserved)
#' all.equal(sort(x[1, ]), sort(x_rand[1, ]))
#' all.equal(sort(x[, 1]), sort(x_rand[, 1]))
#'
#' @export
swapcat <- function(x, n_iter = 1000L, output = c("category", "index"), swaps = "auto", seed = NULL) {
      nullcat(x, method = "swapcat", n_iter = n_iter, output = output, swaps = swaps, seed = seed)
}

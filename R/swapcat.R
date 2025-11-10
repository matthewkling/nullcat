#' Categorical swap (2x2) randomization
#'
#' Categorical generalization of vegan's \code{swap} algorithm.
#'
#' @param x A matrix of categorical data, encoded as integers.
#'   Values should represent category or stratum membership for each cell.
#' @param n_iter Number of swap attempts. Default is 1000.
#' @param output "category" (default) returns randomized matrix;
#'   "index" returns an index matrix describing where original entries moved.
#' @export
swapcat <- function(x, n_iter = 1000L, output = c("category", "index")) {
      output <- match.arg(output)
      x <- as.matrix(x)
      storage.mode(x) <- "integer"
      swapcat_cpp(x, as.integer(n_iter), output)
}

#' Categorical Curveball Randomization
#'
#' Performs a single categorical \emph{curveball} randomization of a matrix,
#' preserving the marginal category multisets of each row and column.
#' This is the categorical generalization of the binary curveball algorithm
#' (Strona et al. 2014), which swaps entries between randomly chosen row pairs
#' while maintaining their category counts.
#'
#' The \code{curvecat()} algorithm treats each cell of \code{x} as belonging
#' to a categorical stratum (integer or factor level). At each iteration, two
#' rows are selected at random and their differing entries are exchanged in
#' random order, ensuring that each row retains its original set of categories.
#' When repeated many times, this procedure generates randomized matrices that
#' preserve row and column category frequencies while randomizing their joint
#' associations.
#'
#' @param x A matrix of categorical data, encoded as integers or factors.
#'   Values should represent category or stratum membership for each cell.
#' @param n_iter Integer specifying the number of randomization iterations
#'   (row-pair trades) to perform. Larger values yield more thorough mixing.
#'   The default is \code{1000}.
#'
#' @details
#' This function provides a convenient R wrapper around the underlying
#' C++ implementation \code{curvecat_cpp()}, which performs the randomization
#' efficiently in compiled code.
#'
#' @return A matrix of the same dimensions and type as \code{x}, with its
#' categorical entries randomized while preserving the marginal category
#' multisets of each row and column.
#'
#' @references
#' Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz, J. (2014).
#' A fast and unbiased procedure to randomize ecological binary matrices
#' with fixed row and column totals. \emph{Nature Communications}, 5, 4114.
#'
#' @seealso
#' \code{\link{trace_curvecat}} for mixing diagnostics,
#' and \code{\link{quantize}} for applying categorical randomization to
#' quantitative community matrices.
#'
#' @examples
#' set.seed(1)
#' m <- matrix(sample(1:4, 20 * 30, replace = TRUE), nrow = 20)
#' m_rand <- curvecat(m, n_iter = 1000)
#' mean(m_rand == m)  # fraction of unchanged cells
#'
#' @export
curvecat <- function(x, n_iter = 1000L) {
      x <- as.matrix(x)
      storage.mode(x) <- "integer"
      curvecat_cpp(x, as.integer(n_iter))
}

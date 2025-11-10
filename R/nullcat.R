
#' Categorical matrix randomization
#'
#' Categorical generalizations of binary community null model algorithms.
#'
#' @param x A matrix of categorical data, encoded as integers.
#'   Values should represent category or stratum membership for each cell.
#' @param method Character specifying the randomization algorithm to use.
#'   Options include:
#'   \itemize{
#'     \item \code{"curvecat"}
#'     \item \code{"swapcat"}
#'   }
#' @param n_iter Number of iterations. Default is 1000.
#' @param output "category" (default) returns randomized matrix;
#'   "index" returns an index matrix describing where original entries moved.
#'
#' @return A matrix of the same dimensions as \code{x}, either randomized
#'   categorical values (\code{output = "category"}) or an integer index
#'   matrix describing the permutation of entries (\code{output = "index"}).
#'
#' @export
nullcat <- function(x,
                    method = c("curvecat", "swapcat"),
                    n_iter = 1000L,
                    output = c("category", "index")) {

      method <- match.arg(method)
      output <- match.arg(output)
      n_iter <- as.integer(n_iter)

      x <- as.matrix(x)
      storage.mode(x) <- "integer"

      switch(
            method,
            curvecat = curvecat(x, n_iter = n_iter, output = output),
            swapcat  = swapcat(x, n_iter = n_iter, output = output)
      )
}


nullcat_methods <- function() c("curvecat", "swapcat")

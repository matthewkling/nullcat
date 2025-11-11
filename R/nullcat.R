
#' Categorical matrix randomization
#'
#' Categorical generalizations of binary community null model algorithms.
#'
#' @param x A matrix of categorical data, encoded as integers.
#'   Values should represent category or stratum membership for each cell.
#' @param method Character specifying the randomization algorithm to use.
#'   Options include:
#'   \itemize{
#'     \item \code{"curvecat"}: categorical analog to `curveball`;
#'          see \link{curvecat} for details.
#'     \item \code{"swapcat"}: categorical analog to `swap`;
#'          see \link{swapcat} for details.
#'     \item \code{"tswapcat"}: categorical analog to `tswap`;
#'          see \link{tswapcat} for details.
#'     \item \code{"r0cat"}: categorical analog to `r0`;
#'          see \link{r0cat} for details.
#'     \item \code{"c0cat"}: categorical analog to `c0`;
#'          see \link{c0cat} for details.
#'   }
#' @param n_iter Number of iterations. Default is 1000.
#'   Ignored for `r0cat` and `col0cat` methods.
#' @param output Character indicating type of result to return:
#'   \itemize{
#'     \item \code{"category"} (default) returns randomized matrix
#'     \item \code{"index"} returns an index matrix describing where original
#'          entries moved.
#'   }
#'
#' @details
#' Note that categorical null models are for cell-level categorical data. Site-level
#' attributes (e.g., land cover) or species-level attributes (e.g., functional
#' traits) should be analyzed using different approaches.
#'
#' @return A matrix of the same dimensions as \code{x}, either randomized
#'   categorical values (when \code{output = "category"}) or an integer index
#'   matrix describing the permutation of entries (when \code{output = "index"}).
#'
#' @export
nullcat <- function(x,
                    method = nullcat_methods(),
                    n_iter = 1000L,
                    output = c("category", "index")) {

      method <- match.arg(method, NULLCAT_METHODS)
      output <- match.arg(output)
      n_iter <- as.integer(n_iter)

      x <- as.matrix(x)
      storage.mode(x) <- "integer"

      switch(
            method,
            "curvecat" = curvecat_cpp(x, n_iter = n_iter, output = output),
            "swapcat"  = swapcat_cpp(x, n_iter = n_iter, output = output),
            "tswapcat"  = tswapcat_cpp(x, n_iter = n_iter, output = output),
            "r0cat"  = r0cat_cpp(x, n_iter = 1, output = output),
            "c0cat"  = c0cat_cpp(x, n_iter = 1, output = output)
      )
}


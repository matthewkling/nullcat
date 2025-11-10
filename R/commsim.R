#' Categorical nullcat commsim (non-sequential)
#'
#' Construct a \code{vegan::commsim} object that uses \code{nullcat()}
#' as a non-sequential null model for categorical / integer matrices.
#' Each simulated matrix is generated independently by applying
#' \code{nullcat()} with \code{n_iter} trades starting from the
#' original matrix.
#'
#' @param n_iter Integer, number of iterations per simulated matrix.
#' @param method Character specifying which nullcat randomization algorithm
#'   to use. See \code{nullcat()} for details.
#' @param output Character, passed to \code{nullcat(output = ...)}.
#'   Typically \code{"category"} (default) or \code{"index"}.
#' @return An object of class \code{"commsim"} suitable for
#'   \code{vegan::nullmodel()}.
#' @export
commsim_cat <- function(n_iter = 1e4,
                        method = "curvecat",
                        output = c("category", "index")) {
      if (!requireNamespace("vegan", quietly = TRUE)) {
            stop("Package 'vegan' is required for commsim_cat(). ",
                 "Please install it with install.packages('vegan').",
                 call. = FALSE)
      }

      output <- match.arg(output)
      n_iter <- as.integer(n_iter)
      if (n_iter < 1L) {
            stop("n_iter must be a positive integer.", call. = FALSE)
      }

      simfun <- function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin, ...) {
            if (!is.matrix(x)) {
                  x <- as.matrix(x)
            }
            if (!is.integer(x)) {
                  storage.mode(x) <- "integer"
            }

            out <- array(NA_integer_, dim = c(nr, nc, n))

            for (k in seq_len(n)) {
                  out[, , k] <- nullcat(x, method = method, n_iter = n_iter, output = output)
            }

            out
      }

      vegan::commsim(
            method = method,
            fun    = simfun,
            binary = FALSE,
            isSeq  = FALSE,
            mode   = "integer"
      )
}




#' Categorical curveball commsim (sequential / Markov chain)
#'
#' Construct a \code{vegan::commsim} object that uses \code{nullcat()}
#' as a *sequential* null model: successive simulated matrices form a
#' Markov chain, and mixing is controlled via the \code{thin} and
#' \code{burnin} arguments to \code{simulate.nullmodel()}.
#'
#' Internally, each simulation "step" advances the chain by
#' \code{thin} curvecat trades.
#'
#' @param method Character specifying which nullcat randomization algorithm
#'   to use. See \code{nullcat()} for details.
#' @param output Character, passed to \code{nullcat(output = ...)}.
#'   Typically \code{"category"} (default) or \code{"index"}.
#' @return An object of class \code{"commsim"} suitable for
#'   \code{vegan::nullmodel()}.
#' @export
commsim_cat_seq <- function(method = "curvecat",
                            output = c("category", "index")) {
      if (!requireNamespace("vegan", quietly = TRUE)) {
            stop("Package 'vegan' is required for commsim_cat_seq(). ",
                 "Please install it with install.packages('vegan').",
                 call. = FALSE)
      }

      output <- match.arg(output)

      simfun <- function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin, ...) {
            thin <- as.integer(thin)
            if (thin < 1L) {
                  stop("thin must be a positive integer for sequential models.", call. = FALSE)
            }

            if (!is.matrix(x)) {
                  x <- as.matrix(x)
            }
            if (!is.integer(x)) {
                  storage.mode(x) <- "integer"
            }

            out   <- array(NA_integer_, dim = c(nr, nc, n))
            state <- x

            for (k in seq_len(n)) {
                  state <- nullcat(state, method = method, n_iter = thin, output = output)
                  out[, , k] <- state
            }

            out
      }

      vegan::commsim(
            method = paste0(method, "_seq"),
            fun    = simfun,
            binary = FALSE,
            isSeq  = TRUE,
            mode   = "integer"
      )
}


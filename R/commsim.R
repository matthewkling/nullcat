#' Categorical curveball commsim (non-sequential)
#'
#' Construct a \code{vegan::commsim} object that uses \code{curvecat()}
#' as a non-sequential null model for categorical / integer matrices.
#' Each simulated matrix is generated independently by applying
#' \code{curvecat()} with \code{n_iter} trades starting from the
#' original matrix.
#'
#' @param n_iter Integer, number of curvecat trades per simulated matrix.
#' @param output Character, passed to \code{curvecat(output = ...)}.
#'   Typically \code{"category"} (default) or \code{"index"}.
#' @param method Character label stored in the \code{commsim} object.
#' @return An object of class \code{"commsim"} suitable for
#'   \code{vegan::nullmodel()}.
#' @export
commsim_curvecat <- function(n_iter = 1e4,
                             output = c("category", "index"),
                             method = "curvecat") {
      if (!requireNamespace("vegan", quietly = TRUE)) {
            stop("Package 'vegan' is required for commsim_curvecat(). ",
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
                  out[, , k] <- curvecat_cpp(x, n_iter, output)
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
#' Construct a \code{vegan::commsim} object that uses \code{curvecat()}
#' as a *sequential* null model: successive simulated matrices form a
#' Markov chain, and mixing is controlled via the \code{thin} and
#' \code{burnin} arguments to \code{simulate.nullmodel()}.
#'
#' Internally, each simulation "step" advances the chain by
#' \code{thin} curvecat trades.
#'
#' @param output Character, passed to \code{curvecat(output = ...)}.
#'   Typically \code{"category"} (default) or \code{"index"}.
#' @param method Character label stored in the \code{commsim} object.
#' @return An object of class \code{"commsim"} suitable for
#'   \code{vegan::nullmodel()}.
#' @export
commsim_curvecat_seq <- function(output = c("category", "index"),
                                 method = "curvecat_seq") {
      if (!requireNamespace("vegan", quietly = TRUE)) {
            stop("Package 'vegan' is required for commsim_curvecat_seq(). ",
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
                  state <- curvecat_cpp(state, thin, output)
                  out[, , k] <- state
            }

            out
      }

      vegan::commsim(
            method = method,
            fun    = simfun,
            binary = FALSE,
            isSeq  = TRUE,
            mode   = "integer"
      )
}


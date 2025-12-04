#' Quantize-based commsim (non-sequential)
#'
#' Construct a [vegan::commsim()] object that uses [quantize()] as a
#' non-sequential null model for numeric community matrices.
#' Each simulated matrix is generated independently by applying
#' [quantize()] with `n_iter` trades (via its internal call to
#' [nullcat()]) starting from the original matrix.
#'
#' @param n_iter Integer, number of iterations (trades) per simulated
#'   matrix. Must be a positive integer. Default is `1e4`.
#' @param ... Arguments passed to [quantize()],
#'   such as `breaks`, `n_strata`, `transform`, `offset`,
#'   `zero_stratum`, `fixed`, `method`, etc. Do **not** supply `x` or
#'   `n_iter` here; these are set internally by `quantize_commsim()`.
#'   See [quantize()] for details.
#'
#' @inheritSection nullcat_commsim Details
#'
#' @seealso [quantize_batch()] if you just want a batch of null matrices
#' without going through **vegan**.
#'
#' @return An object of class `"commsim"` suitable for
#'   [vegan::nullmodel()] and [vegan::oecosimu()].
#'
#' @examplesIf requireNamespace("vegan", quietly = TRUE)
#' \donttest{
#'   library(vegan)
#'
#'   x <- matrix(rexp(50), 10, 5)
#'
#'   cs <- quantize_commsim(
#'     n_strata = 10,
#'     method   = "curvecat",
#'     n_iter   = 1000L
#'   )
#'
#'   nm   <- nullmodel(x, cs)
#'   sims <- simulate(nm, nsim = 999)
#' }
#'
#' @export
quantize_commsim <- function(n_iter = 1e4, ...) {

      require_vegan()
      dots <- validate_dots(list(...))

      n_iter <- as.integer(n_iter)
      if (n_iter < 1L) {
            stop("n_iter must be a positive integer.", call. = FALSE)
      }

      simfun <- function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin, ...) {
            if (!is.matrix(x)) {
                  x <- as.matrix(x)
            }
            if (!is.double(x)) {
                  storage.mode(x) <- "double"
            }

            # one-time prep for this call; pass n_iter through to quantize
            prep <- do.call(
                  quantize_prep,
                  c(list(x = x, n_iter = n_iter), dots)
            )

            out <- array(NA_real_, dim = c(nr, nc, n))

            for (k in seq_len(n)) {
                  out[, , k] <- quantize(prep = prep)
            }

            out
      }

      vegan::commsim(
            method = "quantize",
            fun    = simfun,
            binary = FALSE,
            isSeq  = FALSE,
            mode   = "double"
      )
}



#' Quantile-based quantize commsim (sequential / Markov chain)
#'
#' Construct a [vegan::commsim()] object that uses [quantize()] as a
#' *sequential* null model: successive simulated matrices form a
#' Markov chain in the space of numeric community matrices.
#' Internally, each simulation "step" advances the chain by
#' re-applying [quantize()] to the current matrix using the settings
#' provided via `...`.
#'
#' @param ... Arguments passed to [quantize()], such as `breaks`,
#'   `n_strata`, `transform`, `offset`, `zero_stratum`, `fixed`,
#'   `method`, `n_iter`, etc. Do **not** supply `x` or `n_iter` here;
#'   `x` is provided by vegan and `n_iter` is set internally from
#'   `thin`. See [quantize()] for details.
#'
#' @inheritSection nullcat_commsim_seq Details
#'
#' @return An object of class `"commsim"` suitable for
#'   [vegan::nullmodel()] and [vegan::oecosimu()].
#'
#' @examplesIf requireNamespace("vegan", quietly = TRUE)
#' \donttest{
#'   library(vegan)
#'
#'   x <- matrix(rexp(50), 10, 5)
#'
#'   cs <- quantize_commsim_seq(
#'     n_strata = 5,
#'     method   = "curvecat"
#'   )
#'
#'   nm <- nullmodel(x, cs)
#'
#'   sims <- simulate(
#'     nm,
#'     nsim   = 999,
#'     thin   = 10,    # 10 quantize updates between stored states
#'     burnin = 100    # 100 initial steps discarded
#'   )
#' }
#'
#' @export
quantize_commsim_seq <- function(...) {

      require_vegan()
      dots <- validate_dots(list(...))

      simfun <- function(x, n, nr, nc, rs, cs, rf, cf, s, fill, thin, ...) {
            thin <- as.integer(thin)
            if (thin < 1L) {
                  stop("thin must be a positive integer for sequential models.", call. = FALSE)
            }

            if (!is.matrix(x)) {
                  x <- as.matrix(x)
            }
            storage.mode(x) <- "double"

            out   <- array(NA_real_, dim = c(nr, nc, n))
            state <- x

            for (k in seq_len(n)) {
                  # advance the chain by 'thin' quantize updates
                  for (i in seq_len(thin)) {
                        state <- do.call(quantize, c(list(state), dots))
                  }
                  out[, , k] <- state
            }

            out
      }

      vegan::commsim(
            method = "quantize_seq",
            fun    = simfun,
            binary = FALSE,
            isSeq  = TRUE,
            mode   = "double"
      )
}


validate_dots <- function(dots){
      forbidden <- intersect(names(dots), c("x", "n_iter"))
      if (length(forbidden) > 0L) {
            stop(
                  "The following arguments should not be supplied in `...`: ",
                  paste(forbidden, collapse = ", "),
                  ". They are set internally for this commsim function.",
                  call. = FALSE
            )
      }
      dots
}


#' Nullcat-based commsim (non-sequential)
#'
#' Construct a `vegan::commsim()` object that uses [nullcat()] as a
#' non-sequential null model for categorical / integer matrices.
#' Each simulated matrix is generated independently by applying
#' [nullcat()] with `n_iter` trades starting from the original matrix.
#'
#' @param n_iter Integer, number of iterations (trades) per simulated
#'   matrix. Must be a positive integer. Default is `1e4`.
#' @param method Character specifying which nullcat randomization algorithm
#'   to use. See [nullcat()] and [nullcat_methods()] for details.
#' @param output Character, passed to `nullcat(output = ...)`.
#'   Typically `"category"` (default) or `"index"`.
#'
#' @section Details:
#'
#' This generates a commsim object that is **non-sequential**:
#' each simulated matrix starts from the original matrix and is
#' randomized independently using `n_iter` trades of the chosen
#' `method`.
#'
#' When used via `vegan::simulate.nullmodel()`, the arguments behave as:
#' \itemize{
#'   \item `nsim`: number of simulated matrices to generate.
#'   \item `n_iter` (here, in `nullcat_commsim()`): number of trades per
#'   simulated matrix (controls how strongly each replicate is shuffled).
#'   \item `burnin` and `thin`: are **ignored** for this commsim, because
#'   `isSeq = FALSE` (the simulations are not a Markov chain).
#' }
#'
#' In other words, treat `n_iter` as the tuning parameter for how
#' thoroughly each independent null matrix is randomized.
#'
#' @seealso [nullcat_batch()] if you just want a batch of null matrices
#' without going through **vegan**.
#'
#' @return An object of class `"commsim"` suitable for use with
#'   `vegan::nullmodel()` and `vegan::oecosimu()`.
#'
#' @examples
#' \dontrun{
#'   library(vegan)
#'
#'   x  <- matrix(sample(1:5, 50, replace = TRUE), 10, 5)
#'   cs <- nullcat_commsim(n_iter = 1e4, method = "curvecat")
#'
#'   nm   <- nullmodel(x, cs)
#'   sims <- simulate(nm, nsim = 999)
#' }
#'
#' @export
nullcat_commsim <- function(n_iter = 1e4,
                            method = nullcat_methods(),
                            output = c("category", "index")) {

      require_vegan()
      method <- match.arg(method, NULLCAT_METHODS)
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
                  out[, , k] <- nullcat(
                        x,
                        method = method,
                        n_iter = n_iter,
                        output = output
                  )
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





#' Nullcat-based commsim (sequential / Markov chain)
#'
#' Construct a `vegan::commsim()` object that uses [nullcat()] as a
#' *sequential* null model: successive simulated matrices form a
#' Markov chain. Internally, each simulation "step" advances the chain
#' by `thin` trades of the chosen `method` (e.g. `"curvecat"`), where
#' `thin` is supplied via `vegan::simulate.nullmodel()`` arguments.
#' This is analogous to how sequential swap / curveball null models
#' are used in **vegan**, but extended to categorical data via
#' [nullcat()].
#'
#' @param method Character specifying which nullcat randomization algorithm
#'   to use. See [nullcat()] and [nullcat_methods()] for details.
#' @param output Character, passed to `nullcat(output = ...)`.
#'   Typically `"category"` (default) or `"index"`.
#'
#' @section Details:
#'
#' This model is **sequential**: simulated matrices form a Markov chain. The current matrix is
#' updated in-place by repeated calls to the randomization model, and successive
#' matrices are obtained by advancing the chain.
#'
#' In `vegan::simulate.nullmodel()`, the control arguments behave as:
#' \itemize{
#'   \item `nsim`: number of matrices to *store* from the chain.
#'   \item `thin`: number of trades per step. Each "step" of the chain
#'   applies `thin` trades of the chosen `method` to the current state
#'   before possibly storing it.
#'   \item `burnin`: number of initial steps to perform (each with `thin`
#'   trades) before storing any matrices, i.e. the Markov chain
#'   burn-in.
#' }
#'
#' There is no `n_iter` argument here: mixing is controlled entirely by
#' `thin` (trades per step) and `burnin` (number of initial steps
#' discarded), in the same spirit as sequential swap / curveball models
#' in **vegan**.
#'
#' @return An object of class `"commsim"` suitable for use with
#'   `vegan::nullmodel()` and `vegan::oecosimu()`.
#'
#' @examples
#' \dontrun{
#'   library(vegan)
#'
#'   x  <- matrix(sample(1:5, 50, replace = TRUE), 10, 5)
#'   cs <- nullcat_commsim_seq(method = "curvecat")
#'
#'   nm <- nullmodel(x, cs)
#'
#'   # control the chain with 'thin' and 'burnin'
#'   sims <- simulate(nm, nsim = 999, thin = 100, burnin = 1000)
#' }
#'
#' @export
nullcat_commsim_seq <- function(method = nullcat_methods(),
                                output = c("category", "index")) {

      require_vegan()
      method <- match.arg(method, NULLCAT_METHODS)
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
                  state <- nullcat(
                        state,
                        method = method,
                        n_iter = thin,
                        output = output
                  )
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



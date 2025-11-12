#' Suggest a reasonable n_iter for a randomization
#'
#' Uses trace diagnostics to estimate how many burn-in iterations are
#' needed for a `nullcat` or `quantize` randomization to reach its apparent
#' stationary distribution, given a dataset and randomization method. Uses a
#' "first pre-tail sign-crossing" rule per chain, then returns the maximum
#' across chains. Can be called on a community matrix or a `cat_trace` object.
#'
#' This function uses a “first pre-tail sign-crossing” heuristic to identify burn-in cutoff.
#' This is a simple variant of standard mean-stability tests used in MCMC convergence
#' diagnostics (e.g., Heidelberger & Welch 1983; Geweke 1992; Geyer 1992).
#' It computes the long-run mean based on the "tail window" of the chain, and
#' detects the first iteration at which the trace statistic crosses this
#' long-run mean, indicating that the chain has begun to oscillate around its
#' stationary value. If the chain does not reach the long-run mean before the
#' start of the tail window, the chain is determined not to have reached stationarity,
#' and the function returns \code{NA} with attribute \code{converged = FALSE}.
#'
#' @param trace Either a \code{cat_trace} object (as returned by \code{trace_cat()}), or NULL.
#'   If NULL, arguments to \code{trace_cat()}, including \code{x} and any other relevant
#'   parameters must be supplied via \code{...}
#' @param tail_frac Fraction of the trace (at the end) used as the tail window (default 0.3).
#' @param plot If TRUE, plot the trace, with a vertical line at the suggested value.
#' @param ... Arguments passed to \code{trace_cat()} including  arguments it passes to the
#'   \code{nullcat()} or \code{quantize()} function. Ignored if \code{trace} is non-NULL.
#' @references
#' Heidelberger, P. & Welch, P.D. (1983). Simulation run length control in the presence of an initial transient. Operations Research, 31(6): 1109–1144.
#'
#' Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches to the calculation of posterior moments. In Bayesian Statistics 4, pp. 169–193.
#'
#' Geyer, C.J. (1992). Practical Markov Chain Monte Carlo. Statistical Science, 7(4): 473–483.
#'
#' Feller, W. (1968). An Introduction to Probability Theory and Its Applications, Vol. I. Wiley.
#' @return An integer of class \code{"nullcat_n_iter"} with attributes:
#' \code{n_iter} (numeric or NA), \code{trace} (matrix), \code{steps} (vector),
#' \code{tail_mean} (per-chain), \code{per_chain} (data.frame), \code{converged} (logical).
#'
#' @examples
#' set.seed(1234)
#' x <- matrix(sample(1:5, 2500, replace = T), 50)
#'
#' # call `trace_cat`, then pass result to `suggest_n_iter`:
#' trace <- trace_cat(x = x, fun = "nullcat", n_iter = 1000,
#'                      n_chains = 5, method = "curvecat")
#' suggest_n_iter(trace, tail_frac = 0.3, plot = TRUE)
#'
#' # alternatively, supply `trace_cat` arguments directly to `suggest_n_iter`:
#' x <- matrix(runif(2500), 50)
#' n_iter <- suggest_n_iter(
#'     x = x, n_chains = 5, n_iter = 1000, tail_frac = 0.3,
#'     fun = "quantize", n_strata = 4, fixed = "stratum",
#'     method = "curvecat", plot = T)
#'
#' @export
suggest_n_iter <- function(trace = NULL,
                           tail_frac = 0.3,
                           plot  = FALSE,
                           ...) {

      if(is.null(trace)) trace <- trace_cat(...)

      if (!inherits(trace, "cat_trace")) {
            stop("`trace` must be a cat_trace object (from trace_cat()).")
      }

      tr <- trace$traces
      steps <- trace$steps
      n_steps  <- nrow(tr)
      n_chains <- ncol(tr)

      if (n_steps < 3L || n_chains < 1L) {
            warning("Trace too short to diagnose; increase n_iter or decrease thin.")
            out <- list(n_iter = NA_integer_,
                        trace = tr,
                        steps = steps,
                        tail_mean = numeric(0),
                        per_chain = data.frame(),
                        converged = FALSE)
            class(out) <- "nullcat_n_iter"
            return(out)
      }

      # tail window (indices within rows of `tr`)
      tail_start <- max(1L, n_steps - floor(tail_frac * n_steps) + 1L)

      tail_means      <- numeric(n_chains)
      chain_suggest   <- rep(NA_integer_, n_chains)
      chain_converged <- rep(FALSE,       n_chains)

      for (r in seq_len(n_chains)) {
            chain <- tr[, r]
            mu <- mean(chain[tail_start:n_steps])
            tail_means[r] <- mu

            delta <- chain - mu

            # already exactly at mean at iter 0
            if (is.na(delta[1L]) || delta[1L] == 0) {
                  chain_suggest[r]   <- steps[1L]
                  chain_converged[r] <- TRUE
                  next
            }

            s0 <- sign(delta[1L])

            # must cross BEFORE the tail window starts
            last_pre_tail_idx <- tail_start - 1L
            idx_cross <- NA_integer_

            if (last_pre_tail_idx >= 2L) {
                  for (i in 2L:last_pre_tail_idx) {
                        si <- sign(delta[i])
                        if (!is.na(si) && (si == 0L || si != s0)) {
                              idx_cross <- i
                              break
                        }
                  }
            }

            if (!is.na(idx_cross)) {
                  chain_suggest[r]   <- steps[idx_cross]
                  chain_converged[r] <- TRUE
            } else {
                  chain_suggest[r]   <- NA_integer_
                  chain_converged[r] <- FALSE
            }
      }

      if (!all(chain_converged)) {
            warning(
                  "At least one chain did not cross its tail mean before the tail window.\n",
                  "Increase n_iter (and/or reduce thin) and re-run trace_cat()."
            )
            out <- structure(
                  NA_integer_,
                  trace     = tr,
                  steps     = steps,
                  tail_mean = tail_means,
                  per_chain = data.frame(
                        chain     = seq_len(n_chains),
                        suggested = chain_suggest,
                        converged = chain_converged
                  ),
                  converged = FALSE
            )
            class(out) <- "suggested_n_iter"
            return(out)
      }

      suggested <- max(chain_suggest, na.rm = TRUE)

      if (isTRUE(plot)) {
            plot(trace)
            graphics::abline(v = suggested, lty = 2, col = "black", lwd = 2)
      }

      out <- structure(
            suggested,
            trace     = tr,
            steps     = steps,
            tail_mean = tail_means,
            per_chain = data.frame(
                  chain     = seq_len(n_chains),
                  suggested = chain_suggest,
                  converged = chain_converged
            ),
            converged = TRUE
      )
      class(out) <- "suggested_n_iter"
      out
}


#' @method print suggested_n_iter
#' @export
print.suggested_n_iter <- function(x, ...) {
      cat("suggested_n_iter object\n")
      cat("───────────────────────\n")
      cat("Converged:", attr(x, "converged"), "\n")
      cat("Suggested n iterations:", x, "\n")
      invisible(x)
}

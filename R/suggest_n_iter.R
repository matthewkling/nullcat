#' Suggest a reasonable n_iter for nullcat()
#'
#' Uses trace diagnostics to estimate how many burn-in iterations are
#' needed for the chain to reach its apparent stationary distribution,
#' given your data and randomization method.
#'
#' The diagnostic used is a tail-catching heuristic: burn-in is
#' considered complete when a chosen test statistic first crosses its mean
#' value measured in the final "tail" window of the chain, provided
#' this happens before the tail window. If multiple chains are used
#' (highly recommended), the maximum of burn-in lengths across chains
#' is returned.
#'
#' @inheritParams trace_nullcat
#'
#' @param tail_frac Fraction of the trace (at the end) used to estimate
#'   the stationary mean for each chain. Default is 0.3 (last 30\%).
#'
#' @return
#' A numeric value giving the suggested \code{n_iter} for \code{nullcat()},
#' with attributes \code{trace}, \code{steps}, \code{tail_mean},
#' \code{per_chain}, and \code{converged}. If convergence cannot be
#' diagnosed within \code{n_iter}, the returned value is \code{NA}
#' and \code{converged} is \code{FALSE}.
#'
#' @export
suggest_n_iter <- function(x,
                           method = "curvecat",
                           n_iter = 1000L,
                           thin = NULL,
                           n_chains = 5,
                           stat = "kappa",
                           plot = FALSE,
                           tail_frac = 0.3,
                           n_cores = 1) {

      n_iter <- as.integer(n_iter)

      # Set default thin to give ~100 samples
      if (is.null(thin)) {
            thin <- max(1L, as.integer(n_iter / 100))
      } else {
            thin <- as.integer(thin)
      }

      # Get trace matrix (rows = steps, cols = chains).
      # First row is iteration 0 (x vs x).
      tr <- trace_nullcat(
            x = x,
            method = method,
            n_iter = n_iter,
            thin = thin,
            n_chains = n_chains,
            stat = stat,
            plot = plot,
            n_cores = n_cores
      )

      n_steps  <- nrow(tr)
      n_chains <- ncol(tr)

      if (n_steps == 0L || n_chains == 0L) {
            warning("No trace data; consider increasing n_iter or decreasing thin.")
            return(structure(
                  NA_integer_,
                  trace      = tr,
                  steps      = integer(0),
                  tail_mean  = numeric(0),
                  per_chain  = data.frame(),
                  converged  = FALSE
            ))
      }

      # Iteration numbers corresponding to each trace row:
      # 0, thin, 2*thin, ..., up to ~n_iter
      steps <- seq.int(from = 0L, by = thin, length.out = n_steps)

      tail_means      <- numeric(n_chains)
      chain_suggest   <- rep(NA_integer_, n_chains)
      chain_converged <- rep(FALSE,       n_chains)

      for (r in seq_len(n_chains)) {
            chain <- tr[, r]

            # Estimate stationary value from the tail of this chain
            tail_start <- max(1L, n_steps - floor(tail_frac * n_steps) + 1L)
            tail_vals  <- chain[tail_start:n_steps]
            mu         <- mean(tail_vals)
            tail_means[r] <- mu

            # Deviations from tail mean
            delta <- chain - mu

            # If the initial deviation is exactly zero, treat as already stationary
            if (is.na(delta[1L]) || delta[1L] == 0) {
                  chain_suggest[r]   <- steps[1L]
                  chain_converged[r] <- TRUE
                  next
            }

            # Sign of initial deviation
            s0 <- sign(delta[1L])

            # We only consider sign changes *before* the tail window starts.
            # If the first crossing happens in or after the tail window, we
            # treat that as non-convergence.
            last_pre_tail_idx <- tail_start - 1L

            idx_cross <- NA_integer_

            if (last_pre_tail_idx >= 2L) {
                  for (i in 2L:last_pre_tail_idx) {
                        if (is.na(delta[i])) next
                        si <- sign(delta[i])
                        if (si == 0L || si != s0) {
                              idx_cross <- i
                              break
                        }
                  }
            }

            if (!is.na(idx_cross)) {
                  chain_suggest[r]   <- steps[idx_cross]
                  chain_converged[r] <- TRUE
            } else {
                  # Never crossed the tail mean before the tail window -> non-converged
                  chain_suggest[r]   <- NA_integer_
                  chain_converged[r] <- FALSE
            }
      }

      # If any chain failed to converge, mark overall non-convergence
      if (!all(chain_converged)) {
            warning(
                  "At least one chain did not cross its tail mean before the tail window.\n",
                  "Consider using a larger n_iter (and/or smaller thin)."
            )
            if (plot) {
                  graphics::mtext(
                        "No pre-tail crossing for all chains; increase n_iter",
                        side = 3, col = "red", line = 0.5
                  )
            }

            n <- list(
                  n_iter = NA_integer_,
                  trace = tr,
                  steps = steps,
                  tail_mean = tail_means,
                  per_chain = data.frame(
                        chain = seq_len(n_chains),
                        suggested = chain_suggest,
                        converged = chain_converged
                  ),
                  converged = FALSE)
            class(n) <- "nullcat_n_iter"
            return(n)
      }

      # Overall suggestion: max over per-chain suggestions
      suggested <- max(chain_suggest, na.rm = TRUE)

      if (plot) {
            graphics::abline(v = suggested, lty = 2, col = "black", lwd = 2)
      }

      n <- list(
            n_iter = suggested,
            trace = tr,
            steps = steps,
            tail_mean = tail_means,
            per_chain = data.frame(
                  chain = seq_len(n_chains),
                  suggested = chain_suggest,
                  converged = chain_converged
            ),
            converged = TRUE
      )
      class(n) <- "nullcat_n_iter"
      n
}


#' @method print nullcat_n_iter
#' @export
print.nullcat_n_iter <- function(x, ...) {
      cat("nullcat n_iter object:\n")
      cat("Converged:", x$converged, "\n")
      cat("Suggested n iterations:", x$n_iter, "\n")
      invisible(x)
}

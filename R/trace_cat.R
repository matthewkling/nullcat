
#' Trace diagnostics for categorical randomizations
#'
#' Applies `nullcat()` or `quantize()` to a community matrix, recording
#' a summary statistic at each iteration to help assess mixing on a given dataset.
#'
#' @param x Matrix of categorical data (integers) or quantitative values.
#' @param fun Which function to trace: `"nullcat"` or `"quantize"`.
#' @param n_iter Total number of update iterations to simulate. Default is 1000.
#' @param thin Thinning interval (updates per recorded point). Default ~ `n_iter/100`.
#'   Smaller values increase resolution but increase run time.
#' @param n_chains Number of independent chains to run, to assess consistency (default `5`).
#' @param n_cores Parallel chains (default `1`).
#' @param stat Function that compares `x` to a permuted `x_rand` to quantify their
#'   similarity. Either a function `f(x, x_rand)` returning a scalar, or `NULL`.
#'   If `NULL` (the default), traces use Cohen's kappa for `nullcat()` or Pearson's
#'   correlation for `quantize()`.
#' @param seed Optional integer seed for reproducible traces.
#' @param plot If TRUE, plot the traces.
#' @param ... Arguments to the chosen `fun` (`nullcat()` or `quantize()`),
#'   such as `method`, `n_strata`, `fixed`, etc. If `row_groups` and/or
#'   `col_groups` are supplied (as vectors, or as matrices with one column per
#'   replicate), the traced chains run the grouped randomization exactly as it
#'   is used: every block is advanced by the same number of iterations, and the
#'   statistic is computed on the whole matrix, so the trace (and
#'   [suggest_n_iter()]) gives the number of iterations to run per block. With a
#'   matrix of groups (e.g. tilings), the first column is used. Because every
#'   cell counts equally in the statistic, a few small blocks that mix more
#'   slowly than the rest have little influence on it. Each traced chain costs
#'   as much as one grouped randomization run for `n_iter` iterations (that is,
#'   `n_iter` iterations in every block), so with many blocks it pays not to set
#'   `n_iter` far above the number expected to be needed.
#'
#' @return An object of class `"cat_trace"` with elements:
#' \itemize{
#'   \item `traces`: matrix of size (n_steps+1) x n_chains, including iteration 0
#'   \item `steps`: integer vector of iteration numbers (starting at 0)
#'   \item `fun`, `n_iter`, `thin`, `n_chains`, `n_cores`, `stat_name`, `call`
#'   \item `fun_args`: list of the `...` used (for reproducibility)
#'   \item `n_blocks`: when groups were supplied, the number of blocks in the
#'     traced grouping; otherwise `NULL`
#' }
#' Plotting is available via `plot(cat_trace)`.
#'
#' @examples
#' # nullcat trace
#' set.seed(123)
#' x <- matrix(sample(1:5, 2500, replace = TRUE), 50)
#' tr <- trace_cat(x, n_iter = 1000, n_chains = 5, fun = "nullcat",
#'                 method = "curvecat")
#' plot(tr)
#'
#' # quantize trace
#' x <- matrix(runif(2500), 50)
#' tr <- trace_cat(x, n_iter = 1000, n_chains = 5, fun = "quantize",
#'                 method = "curvecat", n_strata = 3, fixed = "cell")
#' plot(tr)
#' @export
trace_cat <- function(x,
                      fun = c("nullcat","quantize"),
                      n_iter = 1000L,
                      thin = NULL,
                      n_chains = 5L,
                      n_cores = 1L,
                      stat = NULL,
                      seed = NULL,
                      plot = FALSE,
                      ...) {


      fun <- match.arg(fun)

      # choose statistic
      stat_fun <- if (is.null(stat)) {

            if(fun == "nullcat"){
                  stat_name <- "Cohen's kappa"
                  stat <- kappa
            }
            if(fun == "quantize"){
                  stat_name <- "Pearson's r"
                  stat <- function(a, b) cor(as.vector(a), as.vector(b))
            }
            stat
      } else if (is.function(stat)) {
            stat_name <- "Trace statistic"
            stat
      } else {
            stop("if specified, `stat` must be a fuction.")
      }

      core <- trace_chain(
            x0 = x,
            fun = fun,
            n_iter = n_iter,
            thin = thin,
            n_chains = n_chains,
            n_cores = n_cores,
            stat_fun = stat_fun,
            seed = seed,
            ...
      )

      obj <- list(
            traces = core$traces,
            steps = core$steps,
            fun = fun,
            n_iter = as.integer(n_iter),
            thin = if (is.null(thin)) max(1L, as.integer(n_iter/100L)) else as.integer(thin),
            n_chains = as.integer(n_chains),
            n_cores = as.integer(n_cores),
            stat_name = stat_name,
            call = match.call(),
            fun_args = list(...),
            n_blocks = core$n_blocks
      )
      class(obj) <- "cat_trace"

      if (isTRUE(plot)) {
            plot(obj)
      }

      obj
}



# INTERNAL: shared trace core used by trace_cat()
trace_chain <- function(x0,
                        fun = c("nullcat","quantize"),
                        n_iter = 1000L,
                        thin = NULL,
                        n_chains = 5L,
                        n_cores = 1L,
                        stat_fun,
                        seed = NULL,
                        ...) {

      fun <- match.arg(fun)
      n_iter   <- as.integer(n_iter)
      n_chains <- as.integer(n_chains)
      n_cores  <- as.integer(n_cores)

      if (is.null(thin)) thin <- max(1L, as.integer(n_iter / 100L)) else thin <- as.integer(thin)
      steps   <- seq(0L, n_iter, by = thin)[-1L]
      n_steps <- length(steps)
      if (n_steps > 5000L) warning("trace will record ", n_steps, " steps; this may be slow.")

      # With groups, every block is advanced by `thin` iterations between checkpoints
      # (exactly as a grouped randomization runs), and the statistic is computed on
      # the whole matrix. A groups matrix contributes its first column (one tiling).
      dots <- list(...)
      rg <- resolve_groups(dots$row_groups, dots$wt_row, nrow(x0), "row_groups", "wt_row",
                           allow_matrix = TRUE)
      cg <- resolve_groups(dots$col_groups, dots$wt_col, ncol(x0), "col_groups", "wt_col",
                           allow_matrix = TRUE)
      dots[c("row_groups", "col_groups", "wt_row", "wt_col")] <- NULL
      rg <- groups_rep(rg, 1L)
      cg <- groups_rep(cg, 1L)
      grouped <- !is.null(rg) || !is.null(cg)
      n_blocks <- if (grouped) length(block_list(rg, cg, nrow(x0), ncol(x0))) else NULL

      if (fun == "nullcat") {
            update_fun <- function(x) do.call(nullcat, c(list(x, n_iter = thin, output = "category",
                                                              row_groups = rg, col_groups = cg), dots))
            one_chain <- function(i) {
                  x <- x0
                  vals <- numeric(n_steps)
                  for (i in seq_len(n_steps)) {
                        x <- update_fun(x)
                        vals[i] <- stat_fun(x, x0)
                  }
                  vals
            }
      } else { # fun == "quantize" —— optimized path to avoid unneeded overhead computation

            # One-time prep from the initial matrix
            prep  <- do.call(quantize_prep, c(list(x0), dots))
            fixed <- prep$fixed

            # Value pools: global, or one per block when grouped (values never leave their block)
            blocks <- if (grouped) block_list(rg, cg, nrow(x0), ncol(x0)) else NULL
            block_pools <- if (grouped) lapply(blocks, function(b)
                  make_cat_pool(x0[b$i, b$j, drop = FALSE], prep$strata[b$i, b$j, drop = FALSE], fixed)) else NULL
            fill_draw <- function(s_rand) {
                  if (!grouped) {
                        return(fill_from_pool(s = prep$strata, s_rand = s_rand, pool = prep$pool, fixed = fixed))
                  }
                  out <- x0
                  for (k in seq_along(blocks)) {
                        b <- blocks[[k]]
                        out[b$i, b$j] <- fill_from_pool(s = prep$strata[b$i, b$j, drop = FALSE],
                                                        s_rand = s_rand[b$i, b$j, drop = FALSE],
                                                        pool = block_pools[[k]], fixed = fixed)
                  }
                  out
            }
            mode  <- if (fixed == "cell") "index" else "category"

            make_one_chain <- function(i) {
                  # Per-chain state
                  s_cur  <- prep$strata  # current categorical layout
                  x_cur  <- x0 # current quantitative matrix (only used for fixed="cell")

                  # Pre-bind method and thin for speed
                  method <- prep$method
                  nstep  <- thin

                  vals <- numeric(n_steps)
                  for (i in seq_len(n_steps)) {
                        if (mode == "index") {
                              # Advance categories by 'thin' updates and get the composite index permutation
                              idx <- nullcat(s_cur, method = method, n_iter = nstep, output = "index",
                                             row_groups = rg, col_groups = cg)
                              # Move both the categorical layout AND the quantitative values by the same permutation
                              s_cur[] <- s_cur[idx]
                              x_cur[] <- x_cur[idx]
                              # Stat vs. original x0
                              vals[i] <- stat_fun(x_cur, x0)
                        } else {
                              # Advance categorical layout directly
                              s_cur <- nullcat(s_cur, method = method, n_iter = nstep, output = "category",
                                               row_groups = rg, col_groups = cg)
                              # Materialize a quantitative draw using precomputed pool(s)
                              x_draw <- fill_draw(s_cur)
                              vals[i] <- stat_fun(x_draw, x0)
                        }
                  }
                  vals
            }

            one_chain <- make_one_chain
      }

      if (!is.null(seed)) set.seed(as.integer(seed))

      chain_mat <- mc_replicate(n_reps = n_chains, fun = one_chain, n_cores = n_cores)
      chain_mat <- as.matrix(chain_mat)

      out <- matrix(nrow = n_steps + 1L, ncol = n_chains)
      rownames(out) <- paste0("iter", c(0L, steps))
      colnames(out) <- paste0("chain", seq_len(n_chains))
      out[1L, ] <- stat_fun(x0, x0)
      out[2:(n_steps + 1L), ] <- chain_mat

      list(traces = out, steps = c(0L, steps), n_blocks = n_blocks)
}


# Cohen's kappa for categorical data
#
# Computes Cohen's kappa (Cohen, 1960) between two categorical matrices,
# treating each cell as a paired categorical observation.
kappa <- function(a, b) {
      stopifnot(all(dim(a) == dim(b)))
      p_ref <- prop.table(table(as.vector(b)))
      p_chance <- sum(p_ref^2)
      p_obs <- mean(a == b)
      (p_obs - p_chance) / (1 - p_chance)
}




#' @export
plot.cat_trace <- function(x, ...) {
      tr <- x$traces
      steps <- x$steps
      n_chains <- ncol(tr)

      cols <- if (n_chains > 1L) grDevices::rainbow(n_chains) else "black"

      matplot(x = steps, y = tr,
              type = "l", lty = 1, col = cols,
              xlab = "Iteration", ylab = x$stat_name,
              main = paste(x$fun, "mixing traces"),
              ...)

      if (n_chains > 1L) {
            legend("topright", legend = paste("Chain", seq_len(n_chains)),
                   col = cols, lty = 1, bty = "n")
      }

      if(x$stat_name %in% c("Cohen's kappa", "Pearson's r")){
            graphics::abline(h = 0, col = "black", lwd = 1)
      }

      invisible(x)
}


#' @method print cat_trace
#' @export
print.cat_trace <- function(x, digits = 3, ...) {

      if (!inherits(x, "cat_trace")) {
            stop("Object is not of class 'cat_trace'.")
      }

      n_steps <- nrow(x$traces) - 1L  # exclude iteration 0
      n_chains <- ncol(x$traces)
      total_iter <- x$n_iter
      thin <- x$thin
      stat_name <- x$stat_name

      cat("\nCategorical trace diagnostics\n")
      cat("-------------------------------\n")
      cat(sprintf(" Randomization method:   %s\n", x$fun))
      cat(sprintf(" Chains:   %d (%d recorded steps per chain)\n", n_chains, n_steps))
      cat(sprintf(" Iterations: %d total  (thin = %d)\n", total_iter, thin))
      cat(sprintf(" Statistic:  %s\n", stat_name))
      if (!is.null(x$n_blocks)) cat(sprintf(" Grouped:    %d blocks, traced jointly (n_iter per block)\n",
                                            x$n_blocks))
      if (!is.null(x$seed)) cat(sprintf(" Seed:      %d\n", x$seed))
      cat("\n")

      # summarize by chain
      tr <- x$traces
      tail_len <- max(1L, floor(0.2 * nrow(tr)))  # last 20% of trace
      tail_vals <- tr[(nrow(tr) - tail_len + 1L):nrow(tr), , drop = FALSE]
      tail_mean <- apply(tail_vals, 2, mean, na.rm = TRUE)
      tail_sd   <- apply(tail_vals, 2, sd,   na.rm = TRUE)

      df <- data.frame(
            chain = seq_len(n_chains),
            mean = round(tail_mean, digits),
            sd = round(tail_sd, digits)
      )

      print(df, row.names = FALSE)

      cat("\nCall:\n")
      print(x$call)

      invisible(x)
}

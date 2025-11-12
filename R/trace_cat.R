
#' Trace diagnostics for categorical randomizations
#'
#' Applies \code{nullcat()} or \code{quantize()} to a community matrix, recording
#' a summary statistic at each iteration to help assess mixing on a given dataset.
#'
#' @param x Matrix of categorical data (integers) or quantitative values.
#' @param fun Which function to trace: \code{"nullcat"} or \code{"quantize"}.
#' @param n_iter Total number of update iterations to simulate. Default is 1000.
#' @param thin Thinning interval (updates per recorded point). Default ~ \code{n_iter/100}.
#'   Smaller values increase resolution but increase run time.
#' @param n_chains Number of independent chains to run, to assess consistency (default 5).
#' @param n_cores Parallel chains (default 1).
#' @param stat Function that compares \code{x} to a permuted \code{x_rand} to quantify their
#'   similarity. Either a function \code{f(x, x_rand)} returning a scalar, or \code{NULL}.
#'   If \code{NULL} (the default), traces use Cohen's kappa for nullcat or
#'   Pearson's correlation for quantize.
#' @param seed Optional integer seed for reproducible traces.
#' @param plot If TRUE, plot the traces.
#' @param ... Arguments to the chosen \code{fun} (\code{nullcat()} or \code{quantize_null()}),
#'   such as \code{method}, \code{n_strata}, \code{fixed}, etc.
#'
#' @return An object of class \code{"cat_trace"} with elements:
#' \itemize{
#'   \item \code{traces}: matrix of size (n_steps+1) x n_chains, including iteration 0
#'   \item \code{steps}: integer vector of iteration numbers (starting at 0)
#'   \item \code{fun}, \code{n_iter}, \code{thin}, \code{n_chains}, \code{n_cores}, \code{stat_name}, \code{call}
#'   \item \code{fun_args}: list of the \code{...} used (for reproducibility)
#' }
#' Plotting is available via \code{plot(cat_trace)}.
#'
#' @examples
#' # nullcat trace
#' set.seed(123)
#' x <- matrix(sample(1:5, 2500, replace = T), 50)
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
                      fun   = c("nullcat","quantize"),
                      n_iter   = 1000L,
                      thin     = NULL,
                      n_chains = 5L,
                      n_cores  = 1L,
                      stat     = NULL,
                      seed     = NULL,
                      plot     = FALSE,
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
            x0       = x,
            fun   = fun,
            n_iter   = n_iter,
            thin     = thin,
            n_chains = n_chains,
            n_cores  = n_cores,
            stat_fun = stat_fun,
            seed     = seed,
            ...
      )

      obj <- list(
            traces      = core$traces,
            steps       = core$steps,
            fun      = fun,
            n_iter      = as.integer(n_iter),
            thin        = if (is.null(thin)) max(1L, as.integer(n_iter/100L)) else as.integer(thin),
            n_chains    = as.integer(n_chains),
            n_cores     = as.integer(n_cores),
            stat_name   = stat_name,
            call        = match.call(),
            fun_args = list(...)
      )
      class(obj) <- "cat_trace"

      if (isTRUE(plot)) {
            plot(obj)
      }

      obj
}



# INTERNAL: shared trace core used by trace_cat()
trace_chain <- function(x0,
                        fun   = c("nullcat","quantize"),
                        n_iter   = 1000L,
                        thin     = NULL,
                        n_chains = 5L,
                        n_cores  = 1L,
                        stat_fun,
                        seed     = NULL,
                        ...) {

      fun <- match.arg(fun)
      n_iter <- as.integer(n_iter)
      n_chains <- as.integer(n_chains)
      n_cores  <- as.integer(n_cores)

      # default thin ~ n_iter / 100 (≈ 100 samples)
      if (is.null(thin)) {
            thin <- max(1L, as.integer(n_iter / 100L))
      } else {
            thin <- as.integer(thin)
      }

      steps   <- seq(0L, n_iter, by = thin)[-1L]
      n_steps <- length(steps)
      if (n_steps > 5000L) {
            warning("trace will record ", n_steps, " steps; this may be slow.")
      }

      # pick the update kernel for the chosen fun
      update_fun <- switch(
            fun,
            nullcat = function(x) nullcat(x, n_iter = thin, output = "category", ...),
            quantize = function(x) quantize(x, n_iter = thin, ...)
      )

      # one chain (excluding the iteration-0 value)
      one_chain <- function() {
            x <- x0
            vals <- numeric(n_steps)
            for (i in seq_len(n_steps)) {
                  x <- update_fun(x)
                  vals[i] <- stat_fun(x, x0)
            }
            vals
      }

      # reproducible per-chain seeding (uses mc_replicate's internal seeds);
      # setting a seed here ensures mc_replicate’s sampled seeds are reproducible.
      if (!is.null(seed)) {
            set.seed(as.integer(seed))
      }

      chain_mat <- mc_replicate(
            n_reps  = n_chains,
            fun     = one_chain,
            n_cores = n_cores
      )
      chain_mat <- as.matrix(chain_mat)

      # build output with iteration 0 at the top
      out <- matrix(nrow = n_steps + 1L, ncol = n_chains)
      rownames(out) <- paste0("iter", c(0L, steps))
      colnames(out) <- paste0("chain", seq_len(n_chains))

      # row 1: statistic on the original matrix vs itself
      out[1L, ] <- stat_fun(x0, x0)
      out[2:(n_steps + 1L), ] <- chain_mat

      list(traces = out,
           steps  = c(0L, steps))
}





#' Cohen's kappa for categorical data
#'
#' Computes Cohen's kappa between two categorical matrices,
#' treating each cell as a paired categorical observation.
#'
#' @param a,b Matrices or arrays of the same dimensions with categorical entries.
kappa <- function(a, b) {
      stopifnot(all(dim(a) == dim(b)))
      p_ref <- prop.table(table(as.vector(b)))
      p_chance <- sum(p_ref^2)
      p_obs <- mean(a == b)
      (p_obs - p_chance) / (1 - p_chance)
}




cat_mat <- function(nrow = 100, ncol = 50, ncat = 4, prob = NULL){
      matrix(sample(1:ncat, nrow * ncol, replace = T, prob = prob), nrow = nrow)
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
      cat("───────────────────────────────\n")
      cat(sprintf(" Randomization method:   %s\n", x$fun))
      cat(sprintf(" Chains:   %d (%d recorded steps per chain)\n", n_chains, n_steps))
      cat(sprintf(" Iterations: %d total  (thin = %d)\n", total_iter, thin))
      cat(sprintf(" Statistic:  %s\n", stat_name))
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
            mean  = round(tail_mean, digits),
            sd    = round(tail_sd, digits)
      )

      print(df, row.names = FALSE)

      cat("\nCall:\n")
      print(x$call)

      invisible(x)
}

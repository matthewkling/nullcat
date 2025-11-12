

#' Trace diagnostics for nullcat mixing
#'
#' Applies \code{nullcat()} to a categorical matrix, recording a summary statistic
#' at each iteration to help assess mixing on a given dataset.
#'
#' @param x A matrix containing integer or factor values.
#' @param method Character specifying algorithm to use. See \code{nullcat()} for details.
#' @param n_iter Number of successive matrix updates to apply. Default is \code{1000}.
#' @param thin Thinning interval for the test chain. Higher values run faster but
#'   with lower resolution. Default is \code{n_iter / 100}, giving ~100 samples.
#' @param n_chains Number of independent runs to perform, to assess consistency.
#'   Default is \code{5}.
#' @param stat Output trace statistic to compute for each iteration. Either:
#' \itemize{
#'    \item \code{"kappa"} (the default) Cohen's kappa between the original `x` and
#'    the randomized matrix `x_rand`
#'    \item \code{"prop_equal"} proportion of cells in which `x == x_rand`
#'    \item any function taking (`x`, `x_rand`) and returning a numeric scalar
#' }
#' @param plot Logical; if TRUE, plot the trace. Default is FALSE.
#' @param n_cores Maximum number of compute cores to use for parallel processing.
#'   Default is 1. If greater than 1, chains are processed in parallel.
#'
#' @return A matrix of trace stat values, with a row for each reporting iteration
#'    (beginning with iteration 0, the input matrix) and a column for each run.
#'    Row and column names are informative.
#'
#' @export
trace_nullcat <- function(x,
                          method   = "curvecat",
                          n_iter   = 1000L,
                          thin     = NULL,
                          n_chains = 5,
                          stat     = "kappa",
                          plot     = FALSE,
                          n_cores  = 1L) {

      # Choose statistic
      if (stat == "kappa") {
            fun <- kappa
      } else if (stat == "prop_equal") {
            fun <- function(a, b) mean(a == b)
      } else {
            fun <- stat
      }

      # Set default thin to give ~100 samples
      n_iter <- as.integer(n_iter)
      if (is.null(thin)) {
            thin <- max(1L, as.integer(n_iter / 100L))
      } else {
            thin <- as.integer(thin)
      }

      # Reporting steps (excluding iteration 0)
      steps   <- seq(0L, n_iter, by = thin)[-1L]
      n_steps <- length(steps)

      x0 <- x

      # Function to generate one chain's trace (excluding the iteration-0 value)
      one_chain <- function() {
            x    <- x0
            vals <- numeric(n_steps)
            for (i in seq_len(n_steps)) {
                  x <- nullcat(x, method = method, n_iter = thin, output = "category")
                  vals[i] <- fun(x, x0)
            }
            vals
      }

      # Run chains (possibly in parallel)
      # mc_replicate returns an array/matrix with rows = n_steps,
      # cols = n_chains when fun() returns a numeric vector of length n_steps.
      chain_mat <- mc_replicate(
            n_reps  = n_chains,
            fun     = one_chain,
            n_cores = n_cores
      )

      # Ensure it's a matrix with correct orientation
      chain_mat <- as.matrix(chain_mat)
      if (nrow(chain_mat) != n_steps || ncol(chain_mat) != n_chains) {
            stop("Unexpected dimensions from .nullcat_mc_replicate()")
      }

      # Allocate output with an extra row for iteration 0
      out <- matrix(nrow = n_steps + 1L, ncol = n_chains)
      rownames(out) <- paste("iter", c(0L, steps))
      colnames(out) <- paste("chain", seq_len(n_chains))

      # Iteration 0: stat on the original matrix vs itself
      out[1L, ] <- fun(x0, x0)

      # Fill in chain traces (rows 2..n_steps+1)
      out[2:(n_steps + 1L), ] <- chain_mat

      if (plot) {
            cols <- if (n_chains > 1L) grDevices::rainbow(ncol(out)) else "black"
            ylab <- "Trace statistic"
            if (stat == "kappa")      ylab <- "Cohen's kappa"
            if (stat == "prop_equal") ylab <- "Proportion of cells equal"

            matplot(x = c(0L, steps),
                    y = out,
                    type = "l", lty = 1, col = cols,
                    xlab = "Iteration", ylab = ylab,
                    main = paste(method, "mixing traces"))

            if (n_chains > 1L) {
                  legend("topright", legend = paste("Chain", seq_len(ncol(out))),
                         col = cols, lty = 1, bty = "n")
            }
      }

      out
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

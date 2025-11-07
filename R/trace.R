

#' Trace diagnostics for curvecat mixing
#'
#' Applies curvecat to a categorical matrix, recording a summary statistic
#' at each iteration to help assess mixing on a given dataset.
#'
#' @param x A matrix containing integer or factor values.
#' @param max_iter Number of successive curvecat updates to apply. Default is 1000.
#' @param n_runs Number of independent runs to perform. Default is 3.
#' @param stat Output trace statistic to compute for each iteration. Either:
#' \itemize{
#'    \item \code{"kappa"} (the default) Cohen's kappa between the original `x` and
#'    the randomized matrix `x_rand`
#'    \item \code{"prop_equal"} proportion of cells in which `x == x_rand`
#'    \item any function taking (`x`, `x_rand`) and returning a numeric scalar
#' }
#' @param plot Logical; if TRUE, plot the trace. Default is FALSE.
#'
#' @return A matrix of trace stat values, with a row for each iteration
#'    and a column for each run.
#'
#' @export
trace_curvecat <- function(x, max_iter = 1000L, n_runs = 3, stat = "kappa", plot = FALSE) {

      if(stat == "kappa"){
            fun <- kappa
      }else if(stat == "prop_equal"){
            fun <- function(a, b) mean(a == b)
      }else{
            fun <- stat
      }

      x0 <- x
      out <- matrix(nrow = max_iter, ncol = n_runs)
      for(r in 1:n_runs){
            x <- x0
            for (i in seq_len(max_iter)) {
                  x <- curvecat(x, n_iter = 1)
                  out[i,r] <- fun(x, x0)
            }
      }

      if(plot){
            cols <- rainbow(ncol(out))
            ylab <- "Trace statistic"
            if(stat == "kappa") ylab <- "Cohen's kappa"
            if(stat == "prop_equal") ylab <- "Proportion of cells equal"
            matplot(out, type = "l", lty = 1, col = cols,
                    xlab = "Iteration", ylab = ylab,
                    main = "curvecat mixing traces")
            legend("topright", legend = paste("Run", seq_len(ncol(out))),
                   col = cols, lty = 1, bty = "n")
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

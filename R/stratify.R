
#' Bin quantitative data into strata
#'
#' @param x A matrix or vector containing non-negative values.
#' @param breaks Numeric vector of stratum breakpoints.
#' @param n_strata Integer giving the number of strata to split the
#'   data into. Must be 2 or greater. Larger values yield randomizations
#'   with less mixing but higher fidelity to the original marginal
#'   distributions. Default is \code{5}. Ignored unless \code{breaks = NULL}.
#' @param transform A function used to transform the values in
#'   \code{x} before assigning them to \code{n_strata} equal-width
#'   intervals. Examples include \code{sqrt}, \code{log},
#'   \code{rank} (for equal-occupancy strata), etc.;
#'   the default is \code{identity}. If \code{zero_stratum = TRUE}, the
#'   transformation is only applied on nonzero values. The function should
#'   pass NA values. This argument is ignored unless \code{breaks = NULL}.
#' @param offset Numeric value between -1 and 1 (default 0) indicating
#'   how much to shift stratum breakpoints relative to the binwidth (applied
#'   during quantization as: \code{breaks <- breaks + offset * bw}). To
#'   assess sensitivity to stratum boundaries, run \code{quantize()} multiple
#'   times with different offset values. Ignored unless \code{breaks = NULL}.
#' @param zero_stratum Logical indicating whether to segregate zeros into their
#'   own stratum. If \code{FALSE} (the default), zeros will likely be combined
#'   into a stratum that also includes small positive numbers. If \code{breaks} is
#'   specified, zero simply gets added as an additional break; if not, one
#'   of the \code{n_strata} will represent zeros and the others will be nonzero ranges.
#'
#' @return An object the same size as x, with integer values representing
#'   stratum classifications.
#'
#' @examples
#' # Stratify a numeric vector
#' x <- c(0, 0, 0.1, 0.5, 1.2, 3.4, 5.6, 10.2)
#' stratify(x, n_strata = 3)
#'
#' # With transformation
#' stratify(x, n_strata = 3, transform = log1p)
#'
#' # Separate zero stratum
#' stratify(x, n_strata = 3, zero_stratum = TRUE)
#'
#' @export
stratify <- function(x,
                     breaks = NULL,
                     n_strata = 5,
                     transform = identity,
                     offset = 0,
                     zero_stratum = FALSE) {

      if(abs(offset) > 1) stop("`offset` must be between -1 and 1.")
      if(n_strata <= 1 | n_strata != floor(n_strata) | length(n_strata) != 1){
            stop("`n_strata` must be a single integer greater than 1.")
      }

      # Store original dimensions
      x_dim <- dim(x)
      x_dimnames <- dimnames(x)

      # Case 1: User-supplied breaks
      if(!is.null(breaks)){
            if(!length(breaks) > 0 | !is.numeric(breaks)){
                  stop("`breaks` must be a numeric vector with a length of at least 1.")
            }

            if(zero_stratum) {
                  breaks <- c(0, breaks)
            }
            breaks <- unique(sort(c(breaks, -Inf, Inf)))
            s <- as.integer(cut(x, breaks))

            # Restore dimensions
            if(!is.null(x_dim)) {
                  dim(s) <- x_dim
                  dimnames(s) <- x_dimnames
            }

            attr(s, "breaks") <- breaks
            return(s)
      }

      # Case 2: Automatic breaks
      if(zero_stratum){
            # Separate zeros into their own stratum (will be stratum 1)
            zero_mask <- x == 0
            x_nonzero <- x[!zero_mask]

            if(length(x_nonzero) == 0) {
                  # All zeros
                  s <- rep(1L, length(x))
                  if(!is.null(x_dim)) {
                        dim(s) <- x_dim
                        dimnames(s) <- x_dimnames
                  }
                  attr(s, "breaks") <- c(-Inf, 0, Inf)
                  return(s)
            }

            # Check for constant non-zero values
            if(diff(range(x_nonzero, na.rm = TRUE)) == 0) {
                  # All non-zeros are the same value
                  s <- rep(1L, length(x))
                  s[!zero_mask] <- 2L
                  if(!is.null(x_dim)) {
                        dim(s) <- x_dim
                        dimnames(s) <- x_dimnames
                  }
                  attr(s, "breaks") <- c(-Inf, 0, unique(x_nonzero)[1], Inf)
                  return(s)
            }

            # Create n_strata-1 bins for non-zero values
            s_nonzero <- transform(x_nonzero)
            trans <- any(s_nonzero != x_nonzero)

            n_nonzero_strata <- n_strata - 1
            bw <- diff(range(s_nonzero, na.rm = TRUE)) / n_nonzero_strata
            breaks_nonzero <- seq(min(s_nonzero, na.rm = TRUE),
                                  max(s_nonzero, na.rm = TRUE),
                                  length.out = n_nonzero_strata + 1)
            breaks_nonzero <- breaks_nonzero + offset * bw
            breaks_nonzero <- c(-Inf, breaks_nonzero[2:n_nonzero_strata], Inf)

            # Assign strata to non-zero values (will be 1, 2, 3, ...)
            s_nonzero_int <- as.integer(cut(s_nonzero, breaks_nonzero))

            # Initialize full result
            s <- rep(1L, length(x))  # Start with all 1s
            s[!zero_mask] <- as.integer(s_nonzero_int) + 1L  # non-zeros get strata 2, 3, ..., n_strata

            # Restore dimensions
            if(!is.null(x_dim)) {
                  dim(s) <- x_dim
                  dimnames(s) <- x_dimnames
            }

            # Compute effective breaks in original space
            if(trans){
                  # For each stratum 2:n_strata, find max original value
                  bx <- sapply(2:n_strata, function(i) max(x[s == i], na.rm = TRUE))
                  breaks <- c(-Inf, 0, bx[-length(bx)], Inf)
            } else {
                  breaks <- c(-Inf, 0, breaks_nonzero[2:(length(breaks_nonzero)-1)], Inf)
            }

      } else {
            # No zero stratum
            s <- transform(x)
            trans <- any(s != x, na.rm = TRUE)

            # Check for constant values
            if(diff(range(s, na.rm = TRUE)) == 0) {
                  # All values are the same - put everything in stratum 1
                  s[] <- 1L
                  if(!is.null(x_dim)) {
                        dim(s) <- x_dim
                        dimnames(s) <- x_dimnames
                  }
                  attr(s, "breaks") <- c(-Inf, unique(x)[1], Inf)
                  return(s)
            }

            bw <- diff(range(s, na.rm = TRUE)) / n_strata
            breaks <- seq(min(s, na.rm = TRUE),
                          max(s, na.rm = TRUE),
                          length.out = n_strata + 1)
            breaks <- breaks + offset * bw
            breaks <- c(-Inf, breaks[2:n_strata], Inf)

            s[] <- as.integer(cut(s, breaks))

            # Effective breaks in raw space
            if(trans){
                  bx <- sapply(1:n_strata, function(i) max(x[s == i], na.rm = TRUE))
                  breaks <- c(-Inf, bx[1:(n_strata-1)], Inf)
            }
      }

      attr(s, "breaks") <- breaks
      return(s)
}


#' Bin quantitative data into strata
#'
#' @param x A matrix or vector containing non-negative values.
#' @param breaks Numeric vector of stratum breakpoints.
#' @param n_strata Integer giving the number of strata to split the
#'       data into. Must be 2 or greater. Larger values yield randomizations
#'       with less mixing but higher fidelity to the original marginal
#'       distributions. Default is \code{5}. Ignored unless \code{breaks = NULL}.
#' @param transform A function used to transform the values in
#'       \code{x} before assigning them to \code{n_strata} equal-width
#'       intervals. Examples include \code{sqrt}, \code{log}, \code{rank}, etc.;
#'       the default is \code{identity}. If \code{zero_stratum = TRUE}, the
#'       transformation is only applied on nonzero values. The function should
#'       pass NA values. This argument is ignored unless \code{breaks = NULL}.
#' @param offset Numeric value between -1 and 1 (default 0) indicating
#'       how much to shift stratum breakpoints relative to the binwidth (applied
#'       during quantization as: \code{breaks <- breaks + offset * bw}). To
#'       assess sensitivity to stratum boundaries, run \code{quantize()} multiple
#'       times with different offset values. Ignored unless \code{breaks = NULL}.
#' @param zero_stratum Logical indicating whether to segregate zeros into their
#'   own stratum. If \code{FALSE} (the default), zeros will likely be combined
#'   into a stratum that also includes small positive numbers. If \code{breaks} is
#'   specified, zero simply gets added as an additional break; if not, one
#'   of the \code{n_strata} will represent zeros and the others will be nonzero ranges.
#'
#' @return An object the same size as x, with integer values representing
#'   stratum classifications.
#'
#' @export
stratify <- function(x,
                     breaks = NULL,
                     n_strata = 5,
                     transform = identity,
                     offset = 0,
                     zero_stratum = FALSE) {

      if(abs(offset) > 1) stop("`offset` must be between -1 and 1.")

      if(!is.null(breaks)){
            if(zero_stratum) breaks <- c(0, breaks)
            breaks <- unique(sort(c(breaks, -Inf, Inf)))
            x[] <- as.integer(cut(x, breaks))
            attr(x, "breaks") <- breaks
            return(x)
      }

      if(is.null(breaks)){
            if(zero_stratum){
                  n_strata <- n_strata - 1
                  zero <- x == 0
                  x[zero] <- NA
            }

            s <- transform(x)
            bw <- diff(range(s, na.rm = TRUE)) / n_strata
            breaks <- c(-Inf, seq(min(s, na.rm = TRUE) + bw,
                                  max(s, na.rm = TRUE) - bw, bw), Inf)
            breaks <- breaks + offset * bw
            s[] <- as.integer(cut(s, breaks))

            if(zero_stratum){
                  s[zero] <- 0
                  s <- s + 1
            }

            attr(s, "breaks") <- breaks
            return(s)
      }
}

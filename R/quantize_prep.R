#' Prepare stratified null model overhead for quantize()
#'
#' \code{quantize_prep()} precomputes all of the stratification and bookkeeping
#' needed by \code{\link{quantize}()} for a given quantitative community
#' matrix. This is useful when you want to generate many randomizations of the
#' same dataset: the expensive steps (strata assignment, value pools, and
#' arguments for the underlying null model) are computed once, and the resulting
#' object can be passed to \code{quantize(prep = ...)} for fast repeated draws.
#'
#' Internally, \code{quantize_prep()}:
#' \itemize{
#'   \item transforms and stratifies \code{x} into \code{n_strata} numeric
#'     intervals (via \code{\link{stratify}()}),
#'   \item constructs the appropriate value pools given \code{fixed}
#'     (either for the categorical \code{"curvecat"} backend or for binary
#'     vegan methods), and
#'   \item assembles arguments for the underlying null model call
#'     (\code{\link{curvecat}} or \code{\link[vegan]{simulate.nullmodel}}).
#' }
#' The returned object can be reused across calls to \code{\link{quantize}()},
#' \code{\link{quantize_null}()}, or other helpers that accept a \code{prep}
#' argument.
#'
#' @inheritParams quantize
#'
#' @param x Community matrix with sites in rows, species in columns, and
#'   nonnegative quantitative values in cells.
#'   This is the dataset for which stratification and null model overhead
#'   should be prepared.
#'
#' @return A list with class \code{"quantize_prep"} (if you want to set it)
#'   containing the components needed by \code{\link{quantize}()}:
#'   \itemize{
#'     \item \code{strata}: integer matrix of the same dimension as \code{x},
#'       giving the stratum index (1, \dots, \code{n_strata}) for each cell.
#'     \item \code{stratarray}: for binary methods, a 3D array of dimension
#'       \code{n_strata} × \code{nrow(x)} × \code{ncol(x)} containing the
#'       binary incidence matrix for each stratum; \code{NULL} for
#'       \code{method = "curvecat"}.
#'     \item \code{pool}: data structure encoding the quantitative value pools
#'       used during reassignment. For \code{"curvecat"}, this is a list of
#'       per-stratum or per-row/column pools depending on \code{fixed}; for
#'       binary methods, it is a matrix of pre-shuffled values.
#'     \item \code{method}: the null model method used (as in the
#'       \code{method} argument).
#'     \item \code{n_strata}, \code{transform}, \code{offset}, \code{fixed}:
#'       the stratification and reassignment settings used to construct
#'       \code{strata} and \code{pool}.
#'     \item \code{sim_args}: named list of arguments to be passed on to
#'       \code{\link[vegan]{simulate.nullmodel}} (for binary methods) or used
#'       internally by \code{curvecat} (e.g. \code{n_iter}).
#'   }
#'
#'   This object is intended to be passed unchanged to \code{\link{quantize}()}
#'   via its \code{prep} argument.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' comm <- matrix(rexp(50 * 40), nrow = 50,
#'                dimnames = list(paste0("site", 1:50),
#'                                paste0("sp",   1:40)))
#'
#' # prepare overhead for a curvecat-backed stratified null model
#' prep <- quantize_prep(comm, method   = "curvecat",
#'                       n_strata = 5,
#'                       fixed = "row",
#'                       n_iter   = 2000)
#'
#' # fast repeated randomizations using the same prep
#' rand1 <- quantize(prep = prep)
#' rand2 <- quantize(prep = prep)
#'
#' # use a binary vegan method on each stratum
#' prep_bin <- quantize_prep(comm, method = "swap",
#'                           n_strata = 4,
#'                           burnin   = 10000)
#' rand3 <- quantize(prep = prep_bin)
#' }
#'
#' @export
quantize_prep <- function(x,
                          method = nullcat_methods(),
                          fixed = c("stratum", "cell", "row", "col"),
                          breaks = NULL,
                          n_strata  = 5,
                          transform = identity,
                          offset = 0,
                          zero_stratum = FALSE,
                          ...) {

      method <- match.arg(method, NULLCAT_METHODS)
      fixed <- match.arg(fixed)

      # method args
      dots <- list(...)
      sim_args <- dots[! names(dots) %in% names(args)]
      dfts <- list(seed = NULL, burnin = 10000) # defaults
      reqs <- list(nsim = 1, thin = 1) # hard requirements
      sim_args <- c(sim_args, dfts[! names(dfts) %in% names(sim_args)])
      sim_args <- c(reqs, sim_args[! names(sim_args) %in% names(reqs)])

      # [placeholder for error handling for invalid method-fixed combinations]

      # convert to strata
      strata <- stratify(x,
                         breaks = breaks,
                         n_strata = n_strata,
                         transform = transform,
                         offset = offset,
                         zero_stratum = zero_stratum)

      pool <- make_cat_pool(x, strata, fixed = fixed)

      prep <- list(
            strata = strata,
            pool = pool,
            method = method,
            breaks = breaks,
            n_strata = n_strata,
            transform = transform,
            offset  = offset,
            zero_stratum = zero_stratum,
            fixed = fixed,
            sim_args = sim_args
      )
      class(prep) <- c("quantize_prep", "list")
      prep
}


#' @method print quantize_prep
#' @export
print.quantize_prep <- function(x, ...) {
      cat("Quantize preparation object\n")
      cat("Method:", x$method, "\n")
      cat("Strata:", x$n_strata, "\n")
      cat("Priority:", x$priority, "\n")
      cat("Transform:", deparse(x$transform)[1], "\n")
      if (!is.null(x$sim_args$n_iter)) {
            cat("Iterations:", x$sim_args$n_iter, "\n")
      }
      invisible(x)
}


make_cat_pool <- function(x, s, fixed = "stratum") {
      n_strata <- max(s, na.rm = TRUE)

      if (fixed == "cell") {
            pools <- x
      }

      if (fixed == "stratum") {
            pools <- vector("list", n_strata)
            for (k in seq_len(n_strata)) {
                  idx <- which(s == k)
                  pools[[k]] <- x[idx]
            }
      }

      if (fixed == "row") {
            nr <- nrow(x)
            pools <- vector("list", nr)
            for (i in seq_len(nr)) {
                  pools[[i]] <- split(x[i, ], s[i, ])  # list per row, keyed by stratum
            }
      }

      if (fixed == "col") {
            nc <- ncol(x)
            pools <- vector("list", nc)
            for (j in seq_len(nc)) {
                  pools[[j]] <- split(x[, j], s[, j])  # list per col, keyed by stratum
            }
      }

      return(pools)
}

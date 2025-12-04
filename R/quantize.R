
#' Stratified randomization of a quantitative community matrix
#'
#' `quantize()` is a community null model for quantitative community data
#' (e.g. abundance, biomass, or occurrence probability). It works by converting
#' quantitative values into discrete strata, randomizing the stratified matrix
#' using a categorical null model, and reassigning quantitative values within
#' strata according to a specified constraint.
#'
#' This approach provides a framework for preserving row and/or column value
#' distributions in continuous data. When using `fixed = "row"` or
#' `fixed = "col"`, one dimension's value multisets are preserved exactly
#' while the other is preserved at the resolution of strata, approximating a
#' fixed-fixed null model for quantitative data. The number of strata controls
#' the tradeoff between preservation fidelity and randomization strength.
#'
#' By default, `quantize()` will compute all necessary overhead for a
#' given dataset (strata, pools, etc.) internally. For repeated randomization
#' of the same matrix (e.g. to build a null distribution), this overhead can be
#' computed once using [quantize_prep()] and reused by supplying the
#' resulting object via the `prep` argument.
#'
#' @param x Community matrix with sites in rows, species in columns, and
#'   nonnegative quantitative values in cells. Ignored if `prep` is
#'   supplied.
#' @param prep Optional precomputed object returned by
#'   `quantize_prep()`. If supplied, `x` is ignored and all
#'   overhead (stratification, pools, etc.) is taken from `prep`, which
#'   is typically much faster when generating many randomizations of the same
#'   dataset.
#' @param method Character string specifying the null model algorithm.
#'   The default `"curvecat"` uses the categorical curveball algorithm.
#'   See [nullcat()] for alternative options.
#' @param fixed Character string specifying the level at which quantitative
#'   values are held fixed during randomization. One of:
#'   \itemize{
#'     \item `"cell"` (the default; only available when `method = "curvecat"`):
#'       values remain attached to their original cells and move with them during
#'       the categorical randomization. Row and column value
#'       distributions are not preserved, but the mapping between each original
#'       cell and its randomized destination is fixed.
#'     \item `"stratum"`: values are shuffled globally within each stratum,
#'       holding only the overall stratum-level value distribution fixed.
#'     \item `"row"`: values are shuffled within strata separately for each
#'       row, holding each row’s value multiset fixed. Not compatible with all
#'       `method`s.
#'     \item `"col"`: values are shuffled within strata separately for each
#'       column, holding each column’s value multiset fixed.
#' }
#' Note that this interacts with `method`: different null models
#' fix different margins in the underlying binary representation.
#'
#' @inheritParams stratify
#' @inheritParams nullcat
#'
#' @seealso [quantize_batch()] for efficient generation of multiple randomized
#'   matrices; [quantize_commsim()] for integration with `vegan`.
#'
#' @return A randomized version of `x`, with the same dimensions and
#'   dimnames. For `method = "curvecat"`, the quantitative values are
#'   reassigned within strata while preserving row and column stratum
#'   multisets. For binary methods, the result corresponds to applying the
#'   chosen binary null model to each stratum and recombining.
#'
#' @examples
#' # toy quantitative community matrix
#' set.seed(1)
#' comm <- matrix(rexp(50 * 40), nrow = 50,
#'                dimnames = list(paste0("site", 1:50),
#'                                paste0("sp", 1:40)))
#'
#' # default: curvecat-backed stratified randomization
#' rand1 <- quantize(comm)
#'
#' # change stratification and preservation mode
#' rand2 <- quantize(comm, n_strata = 4,
#'                   transform = sqrt,
#'                   fixed  = "row",
#'                   n_iter    = 2000)
#'
#' # use a different randomization algorithm
#' rand3 <- quantize(comm, method = "swapcat", n_iter = 10000)
#'
#' # precompute overhead and reuse for many randomizations
#' prep  <- quantize_prep(comm, method = "curvecat",
#'                        n_strata = 5, fixed = "row")
#' rand4 <- quantize(prep = prep)
#' rand5 <- quantize(prep = prep)
#'
#' @export
#' @rdname quantize
quantize <- function(x = NULL,
                     prep = NULL,
                     method = nullcat_methods(),
                     fixed = c("cell", "stratum", "row", "col"),
                     breaks = NULL,
                     n_strata  = 5,
                     transform = identity,
                     offset = 0,
                     zero_stratum = FALSE,
                     n_iter = 1000,
                     seed = NULL) {

      fixed <- match.arg(fixed)

      if (is.null(prep)) {
            if(is.null(x)) stop("`x` and `prep` cannot both be NULL")
            prep <- quantize_prep(x, method = method, fixed = fixed,
                                  breaks = breaks, n_strata = n_strata, transform = transform,
                                  offset = offset, zero_stratum = zero_stratum, n_iter = n_iter)
      }

      mode <- ifelse(prep$fixed == "cell", "index", "category")
      swaps <- ifelse(prep$fixed == "cell", "alternating", "auto")

      rand_strata <- nullcat(prep$strata, method = prep$method,
                             n_iter = prep$n_iter, output = mode,
                             swaps = swaps, seed = seed)

      with_seed(seed, {
            rand <- fill_from_pool(
                  s = prep$strata,
                  s_rand = rand_strata,
                  pool = prep$pool,
                  fixed = prep$fixed
            )
      })

      return(rand)
}



#' Prepare stratified null model overhead for quantize()
#'
#' `quantize_prep()` precomputes all of the stratification and bookkeeping
#' needed by [quantize()] for a given quantitative community
#' matrix. This is useful when you want to generate many randomizations of the
#' same dataset: the expensive steps (strata assignment, value pools, and
#' arguments for the underlying null model) are computed once, and the resulting
#' object can be passed to `quantize(prep = ...)` for fast repeated draws.
#'
#' Internally, `quantize_prep()`:
#' \itemize{
#'   \item transforms and stratifies `x` into `n_strata` numeric
#'     intervals (via [stratify()]),
#'   \item constructs the appropriate value pools given `fixed`, and
#'   \item assembles arguments for the underlying null model call to [nullcat()].
#' }
#' The returned object can be reused across calls to [quantize()],
#' [quantize_batch()], or other helpers that accept a `prep` argument.
#'
#' @inheritParams quantize
#'
#' @param x Community matrix with sites in rows, species in columns, and
#'   nonnegative quantitative values in cells.
#'   This is the dataset for which stratification and null model overhead
#'   should be prepared.
#'
#' @return A list with class `"quantize_prep"` (if you want to set it)
#'   containing the components needed by `quantize()`:
#'   \itemize{
#'     \item `x`: original quantitative marix `x`,
#'     \item `strata`: integer matrix of the same dimension as `x`,
#'       giving the stratum index (`1:n_strata`) for each cell.
#'     \item `pool`: data structure encoding the quantitative value pools
#'       used during reassignment.
#'     \item `method`: the null model method used (as in the `method` argument).
#'     \item `n_strata`, `transform`, `offset`, `fixed`:
#'       the stratification and reassignment settings used to construct
#'       `strata` and `pool`.
#'     \item `sim_args`: named list of arguments passed to `nullcat()`
#'       (e.g. `n_iter`).
#'   }
#'
#'   This object is intended to be passed unchanged to the `prep` argument of
#'   [quantize()] or [quantize_batch()].
#'
#' @seealso [quantize()], [quantize_batch()]
#'
#' @examples
#' set.seed(1)
#' comm <- matrix(rexp(50 * 40), nrow = 50,
#'                dimnames = list(paste0("site", 1:50),
#'                                paste0("sp", 1:40)))
#'
#' # prepare overhead for a curvecat-backed stratified null model
#' prep <- quantize_prep(comm, method = "curvecat",
#'                       n_strata = 5,
#'                       fixed = "row",
#'                       n_iter = 2000)
#'
#' # fast repeated randomizations using the same prep
#' rand1 <- quantize(prep = prep)
#' rand2 <- quantize(prep = prep)
#'
#' @export
quantize_prep <- function(x,
                          method = nullcat_methods(),
                          fixed = c("cell", "stratum", "row", "col"),
                          breaks = NULL,
                          n_strata  = 5,
                          transform = identity,
                          offset = 0,
                          zero_stratum = FALSE,
                          n_iter = 1000) {

      method <- match.arg(method, NULLCAT_METHODS)
      fixed <- match.arg(fixed)

      # ensure `fixed` is consistent with `method`
      if(method == "r0cat" & fixed == "col" |
         method == "c0cat" & fixed == "row"){
            stop("The selected choices for `fixed` and `method` are incompatible.")
      }

      # convert to strata
      strata <- stratify(x,
                         breaks = breaks,
                         n_strata = n_strata,
                         transform = transform,
                         offset = offset,
                         zero_stratum = zero_stratum)

      pool <- make_cat_pool(x, strata, fixed = fixed)

      prep <- list(
            x = x,
            strata = strata,
            pool = pool,
            method = method,
            breaks = breaks,
            n_strata = n_strata,
            transform = transform,
            offset  = offset,
            zero_stratum = zero_stratum,
            fixed = fixed,
            n_iter = n_iter
      )
      class(prep) <- c("quantize_prep", "list")
      prep
}


#' @method print quantize_prep
#' @export
print.quantize_prep <- function(x, ...) {
      cat("Quantize prep object\n")
      cat("____________________\n")
      cat("Stratification:\n")
      breaks <- signif(attr(x$strata, "breaks"), 2)
      n <- x$n_strata
      strata <- data.frame(stratum = 1:n,
                           from = breaks[1:n],
                           to = breaks[(1:n)+1],
                           freq = as.integer(table(x$strata)))
      print(strata)
      cat("____________________\n")
      cat("Randomization\n")
      cat("   Null algorithm:", x$method, "\n")
      cat("   Quantities fixed:", x$fixed, "\n")

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


fill_from_pool <- function(s, s_rand, pool, fixed = "stratum") {
      out <- matrix(NA_real_, nrow(s), ncol(s))
      dimnames(out) <- dimnames(s)

      fixed <- match.arg(fixed, c("stratum", "row", "col", "cell"))

      if(fixed == "cell"){
            out <- pool
            out[] <- out[s_rand]
            return(out)
      }

      if(fixed == "stratum"){
            n_strata <- length(pool)
            for (k in seq_len(n_strata)) {
                  idx_new <- which(s_rand == k)
                  if (!length(idx_new)) next
                  vals <- pool[[k]]
                  # curvecat preserves global counts per stratum, so lengths should match
                  # just permute without replacement
                  out[idx_new] <- vals[sample.int(length(vals))]
            }
            return(out)
      }

      if(fixed == "row"){
            nr <- nrow(s)
            for (i in seq_len(nr)) {
                  row_pools <- pool[[i]]
                  if (is.null(row_pools)) next

                  # iterate over strata present in this row
                  for (nm in names(row_pools)) {
                        k <- as.integer(nm)
                        idx_new <- which(s_rand[i, ] == k)
                        if (!length(idx_new)) next
                        vals <- row_pools[[nm]]
                        out[i, idx_new] <- vals[sample.int(length(vals))]
                  }
            }
            return(out)
      }

      if(fixed == "col"){
            nc <- ncol(s)
            for (j in seq_len(nc)) {
                  col_pools <- pool[[j]]
                  if (is.null(col_pools)) next

                  for (nm in names(col_pools)) {
                        k <- as.integer(nm)
                        idx_new <- which(s_rand[, j] == k)
                        if (!length(idx_new)) next
                        vals <- col_pools[[nm]]
                        out[idx_new, j] <- vals[sample.int(length(vals))]
                  }
            }
            return(out)
      }
}


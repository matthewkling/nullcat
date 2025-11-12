


#' Stratified randomization of a quantitative community matrix
#'
#' \code{quantize()} is a community null model for quantitative community data
#' (e.g. abundance, biomass, or occurrence probability). It works by converting
#' quantitative values into a small number of categorical strata, randomizing
#' the categorical layout under a chosen categorical null model, and then
#' reassigning quantitative values within each stratum.
#'
#' By default, \code{quantize()} will compute all necessary overhead for a
#' given dataset (strata, pools, etc.) internally. For repeated randomization
#' of the same matrix (e.g. to build a null distribution), this overhead can be
#' computed once using \code{\link{quantize_prep}} and reused by supplying the
#' resulting object via the \code{prep} argument.
#'
#' @param x Community matrix with sites in rows, species in columns, and
#'   nonnegative quantitative values in cells. Ignored if \code{prep} is
#'   supplied.
#' @param prep Optional precomputed object returned by
#'   \code{\link{quantize_prep}}. If supplied, \code{x} is ignored and all
#'   overhead (stratification, pools, etc.) is taken from \code{prep}, which
#'   is typically much faster when generating many randomizations of the same
#'   dataset.
#' @param method Character string specifying the null model algorithm.
#'   The default \code{"curvecat"} uses the categorical curveball algorithm.
#'   See \code{\link{nullcat}} for alternative options.
#' @param fixed Character string specifying the level at which quantitative
#'   values are held fixed during randomization. One of:
#'   \itemize{
#'     \item \code{"cell"} (the default; only available when \code{method = "curvecat"}):
#'       values remain attached to their original cells and move with them during
#'       the categorical randomization. Row and column value
#'       distributions are not preserved, but the mapping between each original
#'       cell and its randomized destination is fixed.
#'     \item \code{"stratum"}: values are shuffled globally within each stratum,
#'       holding only the overall stratum-level value distribution fixed.
#'     \item \code{"row"}: values are shuffled within strata separately for each
#'       row, holding each row’s value multiset fixed.
#'     \item \code{"col"}: values are shuffled within strata separately for each
#'       column, holding each column’s value multiset fixed.
#' }
#' Note that this interacts with \code{method}: different null models
#' fix different margins in the underlying binary representation.
#' @inheritParams stratify
#' @param ... Additional arguments controlling stratification, quantitative
#'   reassignment, and (for vegan methods) the underlying binary null model:
#'   \itemize{
#'     \item \code{n_iter}: For \code{method = "curvecat"}, the number of
#'       categorical curveball iterations (row-pair trades) to perform.
#'       Larger values yield more thorough mixing. If omitted, a data-dependent
#'       default is used.
#'     \item Other arguments are passed on to
#'       \code{\link[vegan]{simulate.nullmodel}} for binary methods, such as
#'       \code{seed} or \code{burnin}. The default \code{burnin} is
#'       \code{10000}. Arguments \code{nsim} and \code{thin} are ignored, as
#'       they are internally set to \code{1}.
#'   }
#'
#' @return A randomized version of \code{x}, with the same dimensions and
#'   dimnames. For \code{method = "curvecat"}, the quantitative values are
#'   reassigned within strata while preserving row and column stratum
#'   multisets. For binary methods, the result corresponds to applying the
#'   chosen binary null model to each stratum and recombining.
#'
#' @examples
#' \donttest{
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
#' # use a vegan binary method on each stratum (here: swap)
#' rand3 <- quantize(comm, method = "swap",
#'                   n_strata = 5,
#'                   burnin   = 10000)
#'
#' # precompute overhead and reuse for many randomizations
#' prep  <- quantize_prep(comm, method = "curvecat",
#'                        n_strata = 5, fixed = "row")
#' rand4 <- quantize(prep = prep)
#' rand5 <- quantize(prep = prep)
#' }
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
                     ...) {

      if (is.null(prep)) {
            if(is.null(x)) stop("`x` and `prep` cannot both be NULL")
            prep <- quantize_prep(x, method = method, fixed = fixed,
                                  breaks = breaks, n_strata = n_strata, transform = transform,
                                  offset = offset, zero_stratum = zero_stratum,
                                  ...)
      }

      n_iter <- prep$sim_args$n_iter
      mode <- ifelse(prep$fixed == "cell", "index", "category")

      rand_strata <- nullcat(prep$strata, method = prep$method, n_iter = n_iter, output = mode)

      rand <- fill_from_pool(
            s = prep$strata,
            s_rand = rand_strata,
            pool = prep$pool,
            fixed = prep$fixed
      )

      return(rand)

}



#' Generate a null distribution using quantize()
#'
#' @param x Community matrix (species × sites, or any numeric matrix).
#' @param n_reps Number of randomizations to generate. Default is `999`.
#' @param stat Optional summary function taking a matrix and returning a numeric
#'        statistic (e.g. `rowSums` with abundance data would give total abundance per site).
#'        If `NULL` (default), the function returns the full set of randomized matrices.
#' @param n_cores Number of compute cores to use for parallel processing. Default is `1`.
#' @param ... Additional arguments passed to `quantize()`
#'        (e.g. `method`, `breaks`, `n_strata`, `transform`, `offset`, `zero_stratum`,
#'        `fixed`, `n_iter`, `burnin`, etc.).
#'
#' @return If stat is NULL: a 3D array (rows × cols × n_reps).
#'   If stat is not NULL: a numeric array of statistic values
#'   (dimensionality will depend on stat).
#'
#' @export
quantize_null <- function(x,
                          n_reps = 999L,
                          stat = NULL,
                          n_cores = 1L,
                          ...) {

      # one-time overhead
      prep <- quantize_prep(as.matrix(x), ...)

      # define per-rep function
      if (is.null(stat)) {
            fun <- function() quantize(prep = prep)
      } else {
            fun <- function() stat(quantize(prep = prep))
      }

      sims <- mc_replicate(n_reps, fun, n_cores = n_cores)

      sims
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


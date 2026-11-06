# Grouped (block) randomization -------------------------------------------------
#
# Row and column groups constrain a randomization so that values are exchanged
# only within blocks: each row-group x column-group submatrix is randomized
# independently with the chosen algorithm. For fixed-margin algorithms this is
# exactly the constrained null: e.g. with row groups, every group x column
# category count is preserved in addition to the usual margins.


# internal -- validate one groups specification. Returns NULL, an integer vector
# of length n, or (if allow_matrix) an integer matrix with n rows, with group
# labels recoded to 1..k within each column.
check_groups <- function(groups, n, arg, allow_matrix = FALSE) {
      if (is.null(groups)) return(NULL)
      dim_name <- if (arg == "row_groups") "nrow(x)" else "ncol(x)"
      if (is.matrix(groups)) {
            if (!allow_matrix) {
                  stop("`", arg, "` must be a vector here. A matrix of groups (one column per ",
                       "replicate) is only accepted by the batch functions.", call. = FALSE)
            }
            if (nrow(groups) != n) {
                  stop("`", arg, "` must have one row per element of ", dim_name, " (", n, ").",
                       call. = FALSE)
            }
            if (anyNA(groups)) stop("`", arg, "` must not contain NA values.", call. = FALSE)
            out <- apply(groups, 2, function(g) match(g, unique(g)))
            dim(out) <- dim(groups)
            storage.mode(out) <- "integer"
            return(out)
      }
      if (!is.atomic(groups) || is.list(groups)) {
            stop("`", arg, "` must be a vector", if (allow_matrix) " or matrix", ".", call. = FALSE)
      }
      if (length(groups) != n) {
            stop("`", arg, "` must have length ", dim_name, " (", n, ").", call. = FALSE)
      }
      if (anyNA(groups)) stop("`", arg, "` must not contain NA values.", call. = FALSE)
      g <- as.vector(groups)
      as.integer(match(g, unique(g)))
}


# internal -- column r of a groups spec (vector specs are returned unchanged).
groups_rep <- function(groups, r) {
      if (is.matrix(groups)) groups[, r] else groups
}


# internal -- check that a groups matrix has one column per replicate.
check_groups_reps <- function(groups, n_reps, arg) {
      if (is.matrix(groups) && ncol(groups) != n_reps) {
            stop("`", arg, "` has ", ncol(groups), " columns but `n_reps` is ", n_reps,
                 "; a groups matrix needs one column per replicate.", call. = FALSE)
      }
      invisible(NULL)
}


# internal -- list of blocks (row indices, column indices) for given groups.
block_list <- function(row_groups, col_groups, nr, nc) {
      rg <- if (is.null(row_groups)) list(seq_len(nr)) else split(seq_len(nr), row_groups)
      cg <- if (is.null(col_groups)) list(seq_len(nc)) else split(seq_len(nc), col_groups)
      blocks <- vector("list", length(rg) * length(cg))
      k <- 0L
      for (i in rg) for (j in cg) {
            k <- k + 1L
            blocks[[k]] <- list(i = i, j = j)
      }
      blocks
}


# internal -- apply `f` (a function of a submatrix returning a same-sized
# matrix) to each block of `x`, reassembling the result. With
# `output = "index"`, `f` returns linear indices into the submatrix, which are
# translated into linear indices into `x`.
apply_blocks <- function(x, row_groups, col_groups, f, output = "category") {
      out <- if (output == "index") matrix(seq_along(x), nrow(x), ncol(x)) else x
      for (b in block_list(row_groups, col_groups, nrow(x), ncol(x))) {
            sub <- x[b$i, b$j, drop = FALSE]
            r <- f(sub)
            if (output == "index") {
                  global <- matrix(seq_along(x), nrow(x), ncol(x))[b$i, b$j, drop = FALSE]
                  out[b$i, b$j] <- global[r]
            } else {
                  out[b$i, b$j] <- r
            }
      }
      out
}


# internal -- turn a (deprecated) pairwise weight matrix into groups.
# Symmetric weighted pair selection only changes how fast a chain mixes, not
# what it converges to; the only lasting constraint is which rows can ever
# exchange values, i.e. the connected components of the nonzero weights.
weights_to_groups <- function(W, n, arg) {
      if (!is.matrix(W) || !is.numeric(W)) stop(arg, " must be a numeric matrix.", call. = FALSE)
      if (nrow(W) != ncol(W)) stop(arg, " must be a square matrix.", call. = FALSE)
      if (nrow(W) != n) {
            dim_name <- if (arg == "wt_row") "nrow" else "ncol"
            stop(arg, " must be ", n, " x ", n, " (matching ", dim_name, "(x)).", call. = FALSE)
      }
      if (anyNA(W)) stop(arg, " must not contain NA values.", call. = FALSE)
      if (any(W < 0)) stop(arg, " must contain only non-negative values.", call. = FALSE)

      A <- W > 0 | t(W) > 0
      diag(A) <- FALSE
      if (!any(A)) {
            stop("`", arg, "` has no positive off-diagonal weights, so no ",
                 if (arg == "wt_row") "rows" else "columns",
                 " could exchange values.", call. = FALSE)
      }
      groups <- connected_components(A)
      k <- max(groups)
      new_arg <- if (arg == "wt_row") "row_groups" else "col_groups"
      if (k == 1L) {
            warning("`", arg, "` is deprecated and has no effect here: its nonzero weights connect all ",
                    if (arg == "wt_row") "rows" else "columns",
                    ", and weighting which pairs exchange values does not change the null ",
                    "distribution a randomization converges to (only how fast it gets there). ",
                    "Use `", new_arg, "` to constrain randomization within groups.", call. = FALSE)
            return(NULL)
      }
      warning("`", arg, "` is deprecated; use `", new_arg, "`. The weights were converted to ", k,
              " groups (the connected components of their nonzero entries); only this ",
              "grouping affects the null distribution, not the weight values.", call. = FALSE)
      groups
}


# internal -- connected components of a symmetric logical adjacency matrix.
connected_components <- function(A) {
      n <- nrow(A)
      comp <- integer(n)
      k <- 0L
      for (s in seq_len(n)) {
            if (comp[s] > 0L) next
            k <- k + 1L
            comp[s] <- k
            frontier <- s
            while (length(frontier)) {
                  nb <- which(colSums(A[frontier, , drop = FALSE]) > 0 & comp == 0L)
                  comp[nb] <- k
                  frontier <- nb
            }
      }
      comp
}


# internal -- reconcile groups and deprecated weights for one margin.
resolve_groups <- function(groups, wt, n, arg_groups, arg_wt, allow_matrix = FALSE) {
      if (!is.null(wt)) {
            if (!is.null(groups)) {
                  stop("Supply only one of `", arg_groups, "` and `", arg_wt, "`.", call. = FALSE)
            }
            return(weights_to_groups(wt, n, arg_wt))
      }
      check_groups(groups, n, arg_groups, allow_matrix = allow_matrix)
}


#' Mean group size for each row
#'
#' For a grouping vector or a matrix of groupings (such as the output of
#' [tile_groups()], with one column per randomization), returns the average
#' number of members in each element's group. With a groups matrix, this is the
#' average over columns. It describes the size of the pool each row (e.g. site)
#' is randomized within, and depends only on the grouping, so it can be checked
#' before running any randomizations, for example while choosing a tile
#' `spacing`. Rows that are often in very small groups are rarely changed by a
#' fixed-margin randomization; a row alone in its group is never changed.
#'
#' @param groups A vector of group labels, or a matrix of group labels with one
#'   column per randomization.
#'
#' @return A numeric vector with one value per element (row) of `groups`.
#'
#' @seealso [tile_groups()]
#'
#' @examples
#' group_sizes(c("a", "a", "b", "c", "c", "c"))
#'
#' set.seed(1)
#' coords <- cbind(runif(200, 0, 10), runif(200, 0, 10))
#' g <- tile_groups(coords, spacing = 3, n = 50, lonlat = FALSE)
#' summary(group_sizes(g))
#'
#' @export
group_sizes <- function(groups) {
      if (is.null(groups)) stop("`groups` must not be NULL.", call. = FALSE)
      if (anyNA(groups)) stop("`groups` must not contain NA values.", call. = FALSE)
      size1 <- function(g) {
            g <- match(g, unique(g))
            as.numeric(tabulate(g)[g])
      }
      if (is.matrix(groups)) {
            sizes <- vapply(seq_len(ncol(groups)), function(r) size1(groups[, r]),
                            numeric(nrow(groups)))
            return(rowMeans(matrix(sizes, nrow(groups))))
      }
      size1(as.vector(groups))
}

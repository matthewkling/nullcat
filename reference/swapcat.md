# Categorical swap randomization (swapcat)

Categorical generalization of the binary 2x2 swap algorithm to matrices
of categorical data. This function is a convenience wrapper around
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
with `method = "swapcat"`.

## Usage

``` r
swapcat(
  x,
  n_iter = 1000L,
  output = c("category", "index"),
  swaps = "auto",
  row_groups = NULL,
  col_groups = NULL,
  wt_row = NULL,
  wt_col = NULL,
  seed = NULL
)
```

## Arguments

- x:

  A matrix of categorical data, encoded as integers. Values should
  represent category or stratum membership for each cell.

- n_iter:

  Number of iterations. Default is 1000. Larger values yield more
  thorough mixing. Ignored for non-sequential methods. Minimum burn-in
  times can be estimated with
  [`suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md).
  When `row_groups` or `col_groups` is used, `n_iter` iterations are run
  separately on every block. The number needed can be more or fewer than
  for the full matrix, so estimate it by passing the same groups to
  [`suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md).

- output:

  Character indicating type of result to return:

  - `"category"` (default) returns randomized matrix

  - `"index"` returns an index matrix describing where original entries
    (a.k.a. "tokens") moved. Useful mainly for testing, and for
    applications like
    [`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md)
    that care about token tracking in addition to generic integer
    categories.

- swaps:

  Character string controlling the direction of token movement. Only
  used when method is `"curvecat"`, `"swapcat"`, or `"tswapcat"`.
  Affects the result only when `output = "index"`, otherwise it only
  affects computation speed. Options include:

  - `"vertical"`: Tokens move between rows (stay within columns).

  - "`horizontal"`: Tokens move between columns (stay within rows).

  - `"alternating"`: Tokens move in both dimensions, alternating between
    vertical and horizontal swaps. Provides full 2D mixing without
    preserving either row or column token sets.

  - `"auto"` (default): For `output = "category"`, automatically selects
    the fastest option based on matrix dimensions. For
    `output = "index"`, defaults to `"alternating"` for full mixing.

- row_groups, col_groups:

  Optional vectors assigning each row (for `row_groups`) or column (for
  `col_groups`) of `x` to a group. Labels can be integers, characters,
  or factor levels. When supplied, the matrix is split into blocks
  (every combination of a row group and a column group) and each block
  is randomized independently with `method`, so values are only ever
  exchanged within a block. For algorithms that fix row and column
  margins, this additionally fixes the category counts of every row
  group within each column (and every column group within each row): for
  example, grouping sites into regions gives a null model in which each
  region keeps its own species pool. Rows or columns that are alone in
  their group cannot be changed by fixed-margin algorithms. To randomize
  within a different grouping in each replicate (e.g. random spatial
  tilings from
  [`tile_groups()`](https://matthewkling.github.io/nullcat/reference/tile_groups.md)),
  pass a groups matrix to
  [`nullcat_batch()`](https://matthewkling.github.io/nullcat/reference/nullcat_batch.md)
  or
  [`quantize_batch()`](https://matthewkling.github.io/nullcat/reference/quantize_batch.md).

- wt_row, wt_col:

  Deprecated. Experimental pairwise weight matrices previously used to
  bias which pairs of rows or columns exchange values. But weighting
  pair selection changes only how quickly a randomization mixes, not the
  null distribution it converges to; the only lasting effect of a weight
  matrix is to prevent exchanges between rows (or columns) that are not
  connected, directly or through intermediaries, by nonzero weights.
  Weight matrices are therefore now converted to `row_groups` /
  `col_groups` given by the connected components of their nonzero
  entries, with a warning. Use `row_groups` / `col_groups` instead.

- seed:

  Integer used to seed random number generator, for reproducibility.

## Value

A matrix of the same dimensions as `x`, either randomized categorical
values (when `output = "category"`) or an integer index matrix
describing the permutation of entries (when `output = "index"`).

## Details

The swapcat algorithm attempts random 2x2 swaps of the form:

    a b        b a
    b a   <->  a b

where \\a\\ and \\b\\ are distinct categories. These swaps preserve the
multiset of categories in each row and column. With only two categories
present, `swapcat()` reduces to the behavior of the standard binary swap
algorithm.

## References

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. *Ecology*, 81(9), 2606–2621.

## See also

[`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md)
for an algorithm that produces equivalent results with better
computational efficiency.

## Examples

``` r
set.seed(123)
x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
x_rand <- swapcat(x, n_iter = 1000)
```

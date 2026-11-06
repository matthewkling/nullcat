# Trial-swap categorical randomization (tswapcat)

The trial-swap ("tswap") algorithm is a fixed-fixed randomization that
repeatedly attempts random 2x2 swaps until a valid one is found in each
iteration, reducing the number of wasted draws compared to the simple
swap. `tswapcat()` extends this logic to categorical matrices.

## Usage

``` r
tswapcat(
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
  separately on every block, so it should be sized for a block rather
  than for the full matrix.

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

Unlike
[`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md)
and
[`swapcat()`](https://matthewkling.github.io/nullcat/reference/swapcat.md),
whose symmetric moves converge to a uniform draw from all matrices with
the observed row and column category multisets, redrawing until a valid
swap is found biases sampling toward matrices with many available swaps
(Miklos & Podani 2004). The stationary distribution of `tswapcat()` is
therefore not uniform, and it is provided mainly for methodological
comparison.

## References

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. *Ecology*, 81(9), 2606-2621.

Miklos, I. & Podani, J. (2004). Randomization of presence-absence
matrices: comments and new algorithms. *Ecology*, 85(1), 86-92.

## See also

[`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md)
for an efficient algorithm that samples matrices uniformly.

## Examples

``` r
set.seed(123)
x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
x_rand <- tswapcat(x, n_iter = 1000)
```

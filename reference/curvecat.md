# Categorical curveball randomization (curvecat)

Categorical generalization of the binary curveball algorithm (Strona et
al. 2014) to matrices of categorical data. This function is a
convenience wrapper around
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
with `method = "curvecat"`.

## Usage

``` r
curvecat(
  x,
  n_iter = 1000L,
  output = "category",
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

The curvecat algorithm randomizes a categorical matrix while keeping the
category multisets of each row and column fixed. In other words, the
permuted matrix has the same set of integer values in every row and
every column as the original matrix, but they are permuted. It operates
on pairs of rows at a time, grouping differing entries by unordered
category pairs and redistributing the orientation of those pairs while
preservingn the multiset of categories within each row. When there are
only two categories, `curvecat()` reduces to the behavior of the
original binary curveball algorithm (Strona et al. 2014) applied to a
0/1 matrix.

## References

Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz,
J. (2014). A fast and unbiased procedure to randomize ecological binary
matrices with fixed row and column totals. *Nature Communications*, 5,
4114.

## See also

[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md),
[`nullcat_methods()`](https://matthewkling.github.io/nullcat/reference/nullcat_methods.md)

## Examples

``` r
# Create a categorical matrix
set.seed(123)
x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

# Randomize preserving row and column category multisets
x_rand <- curvecat(x, n_iter = 1000)

# Verify margins are preserved
all.equal(sort(x[1, ]), sort(x_rand[1, ])) # row multisets preserved
#> [1] TRUE
all.equal(sort(x[, 1]), sort(x_rand[, 1])) # column multisets preserved
#> [1] TRUE
```

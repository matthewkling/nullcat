# Categorical matrix randomization

Randomize binary or categorical community matrices using categorical
generalizations of binary community null model algorithms. Optionally
constrain randomization to occur within groups of rows (e.g. regions)
and/or columns (e.g. clades).

## Usage

``` r
nullcat(
  x,
  method = nullcat_methods(),
  n_iter = 1000L,
  output = c("category", "index"),
  swaps = c("auto", "vertical", "horizontal", "alternating"),
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

- method:

  Character specifying the randomization algorithm to use. Options
  include the following; see details and linked functions for more info.

  - `"curvecat"`: categorical analog to `curveball`; see
    [`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md)
    for details.

  - `"swapcat"`: categorical analog to `swap`; see
    [`swapcat()`](https://matthewkling.github.io/nullcat/reference/swapcat.md)
    for details.

  - `"tswapcat"`: categorical analog to `tswap`; see
    [`tswapcat()`](https://matthewkling.github.io/nullcat/reference/tswapcat.md)
    for details.

  - `"r0cat"`: categorical analog to `r0`; see
    [`r0cat()`](https://matthewkling.github.io/nullcat/reference/r0cat.md)
    for details.

  - `"c0cat"`: categorical analog to `c0`; see
    [`c0cat()`](https://matthewkling.github.io/nullcat/reference/c0cat.md)
    for details.

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

`curvecat`, `swapcat`, and `tswapcat` are sequential algorithms that
hold category multisets fixed in every row and column. `curvecat` and
`swapcat` propose symmetric moves, so with enough iterations they
converge to the same null distribution: a uniform draw from all matrices
with the observed row and column category multisets. `curvecat` is the
most efficient of the two (i.e. fewest steps to become fully mixed).
`tswapcat` repeatedly redraws until it finds a valid swap; this makes
each iteration productive, but biases the sampling toward matrices with
many available swaps, so its stationary distribution is not uniform (cf.
Miklós & Podani 2004). `swapcat` and `tswapcat` are thus useful mainly
for methodological comparison.

The `r0cat` algorithm holds category multisets fixed in rows but not
columns, while `c0cat` does the opposite.

Note that categorical null models are for cell-level categorical data.
Site-level attributes (e.g., land cover) or species-level attributes
(e.g., functional traits) should be analyzed using different approaches.
See vignette for details.

## References

Miklós, I., & Podani, J. (2004). Randomization of presence–absence
matrices: comments and new algorithms. *Ecology*, 85(1), 86–92.

## See also

[`nullcat_batch()`](https://matthewkling.github.io/nullcat/reference/nullcat_batch.md)
for efficient generation of multiple randomized matrices;
[`nullcat_commsim()`](https://matthewkling.github.io/nullcat/reference/nullcat_commsim.md)
for integration with `vegan`;
[`tile_groups()`](https://matthewkling.github.io/nullcat/reference/tile_groups.md)
for random spatial groupings.

## Examples

``` r
# Create a categorical matrix
set.seed(123)
x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

# Randomize using curvecat method (preserves row & column margins)
x_rand <- nullcat(x, method = "curvecat", n_iter = 1000)

# Check that row multisets are preserved
all.equal(sort(x[1, ]), sort(x_rand[1, ]))
#> [1] TRUE

# Randomize within two regions of five sites each
region <- rep(c("north", "south"), each = 5)
x_reg <- nullcat(x, method = "curvecat", n_iter = 1000, row_groups = region)

# Each region keeps its own category counts in every column
all.equal(rowsum((x == 1) * 1, region), rowsum((x_reg == 1) * 1, region))
#> [1] TRUE
```

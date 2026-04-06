# Categorical matrix randomization

Randomize binary or categorical community matrices using categorical
generalizations of binary community null model algorithms. Optionally
constrain mixing using spatial (row) and taxonomic (column) weights.

## Usage

``` r
nullcat(
  x,
  method = nullcat_methods(),
  n_iter = 1000L,
  output = c("category", "index"),
  swaps = c("auto", "vertical", "horizontal", "alternating"),
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
    When `wt_row` or `wt_col` is supplied, defaults to the appropriate
    direction, or `"alternating"` if both are supplied.

- wt_row:

  An optional square numeric matrix of non-negative weights controlling
  which pairs of rows are likely to exchange tokens during
  randomization. Must be `nrow(x)` by `nrow(x)`. This enables spatially
  or trait-constrained null models where nearby or similar sites
  exchange tokens more frequently.

  Values are treated as relative weights (not probabilities) and are
  normalized internally. The diagonal is ignored. The matrix should be
  symmetric. Only supported for sequential methods (`curvecat`,
  `swapcat`, `tswapcat`).

  When both `wt_row` and `wt_col` are supplied, `swaps` is forced to
  `"alternating"`, producing a Gibbs-like sweep that applies each weight
  matrix on its respective margin in alternation.

- wt_col:

  An optional square numeric matrix of non-negative weights controlling
  which pairs of columns are likely to exchange tokens during
  randomization. Must be `ncol(x)` by `ncol(x)`. See `wt_row` for
  details on weight interpretation.

- seed:

  Integer used to seed random number generator, for reproducibility.

## Value

A matrix of the same dimensions as `x`, either randomized categorical
values (when `output = "category"`) or an integer index matrix
describing the permutation of entries (when `output = "index"`).

## Details

`curvecat`, `swapcat`, and `tswapcat` are sequential algorithms that
hold category multisets fixed in every row and column. These three
algorithms typically reach the same stationary distribution. They differ
primarily in efficiency, with `curvecat` being the most efficient (i.e.
fewest steps to become fully mixed); `swapcat` and `tswapcat` are thus
useful mainly for methodological comparison.

The `r0cat` algorithm holds category multisets fixed in rows but not
columns, while `c0cat` does the opposite.

Note that categorical null models are for cell-level categorical data.
Site-level attributes (e.g., land cover) or species-level attributes
(e.g., functional traits) should be analyzed using different approaches.
See vignette for details.

## See also

[`nullcat_batch()`](https://matthewkling.github.io/nullcat/reference/nullcat_batch.md)
for efficient generation of multiple randomized matrices;
[`nullcat_commsim()`](https://matthewkling.github.io/nullcat/reference/nullcat_commsim.md)
for integration with `vegan`.

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

# Spatially constrained randomization using row weights
coords <- cbind(runif(10), runif(10))
d <- as.matrix(dist(coords))
W <- exp(-d / 0.3)  # Gaussian distance decay
x_spatial <- nullcat(x, method = "curvecat", n_iter = 1000, wt_row = W)

# Dual-margin weighting (Gibbs-like alternating)
W_row <- exp(-as.matrix(dist(cbind(runif(10), runif(10)))) / 0.3)
W_col <- exp(-as.matrix(dist(cbind(runif(10), runif(10)))) / 0.3)
x_dual <- nullcat(x, method = "curvecat", n_iter = 1000,
                  wt_row = W_row, wt_col = W_col)
```

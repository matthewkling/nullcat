# Categorical matrix randomization

Categorical generalizations of binary community null model algorithms.

## Usage

``` r
nullcat(
  x,
  method = nullcat_methods(),
  n_iter = 1000L,
  output = c("category", "index"),
  swaps = c("auto", "vertical", "horizontal", "alternating"),
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

  - `"curvecat"`: categorical analog to \`curveball\`; see
    [curvecat](https://matthewkling.github.io/nullcat/reference/curvecat.md)
    for details.

  - `"swapcat"`: categorical analog to \`swap\`; see
    [swapcat](https://matthewkling.github.io/nullcat/reference/swapcat.md)
    for details.

  - `"tswapcat"`: categorical analog to \`tswap\`; see
    [tswapcat](https://matthewkling.github.io/nullcat/reference/tswapcat.md)
    for details.

  - `"r0cat"`: categorical analog to \`r0\`; see
    [r0cat](https://matthewkling.github.io/nullcat/reference/r0cat.md)
    for details.

  - `"c0cat"`: categorical analog to \`c0\`; see
    [c0cat](https://matthewkling.github.io/nullcat/reference/c0cat.md)
    for details.

- n_iter:

  Number of iterations. Default is 1000. Larger values yield more
  thorough mixing. Ignored for non-sequential methods. Minimum burn-in
  times can be estimated with
  [suggest_n_iter](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md).

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
  used when method is \`curvecat\`, \`swapcat\`, or \`tswapcat\`.
  Affects the result only when `output = "index"`, otherwise it only
  affects computation speed. Options include:

  - `"vertical"`: Tokens move between rows (stay within columns).

  - `"horizontal"`: Tokens move between columns (stay within rows).

  - `"alternating"`: Tokens move in both dimensions, alternating between
    vertical and horizontal swaps. Provides full 2D mixing without
    preserving either row or column token sets.

  - `"auto"` (default): For `output = "category"`, automatically selects
    the fastest option based on matrix dimensions. For
    `output = "index"`, defaults to `"alternating"` for full mixing.

- seed:

  Integer used to seed random number generator, for reproducibility.

## Value

A matrix of the same dimensions as `x`, either randomized categorical
values (when `output = "category"`) or an integer index matrix
describing the permutation of entries (when `output = "index"`).

## Details

\`curvecat\`, \`swapcat\`, and \`tswapcat\` are sequential algorithms
that hold category multisets fixed in every row and column. These three
algorithms typically reach the same stationary distribution. They differ
primarily in efficiency, with \`curvecat\` being the most efficient
(i.e. fewest steps to become fully mixed); \`swapcat\` and \`tswapcat\`
are thus useful mainly for methodological comparison.

The \`r0cat\` algorithm holds category multisets fixed in rows but not
columns, while \`c0cat\` does the opposite.

Note that categorical null models are for cell-level categorical data.
Site-level attributes (e.g., land cover) or species-level attributes
(e.g., functional traits) should be analyzed using different approaches.

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

# Get index output showing where each cell moved
idx <- nullcat(x, method = "curvecat", n_iter = 1000, output = "index")

# Use different methods
x_swap <- nullcat(x, method = "swapcat", n_iter = 1000)
x_r0 <- nullcat(x, method = "r0cat")

# Use with a seed for reproducibility
x_rand1 <- nullcat(x, method = "curvecat", n_iter = 1000, seed = 42)
x_rand2 <- nullcat(x, method = "curvecat", n_iter = 1000, seed = 42)
identical(x_rand1, x_rand2)
#> [1] TRUE
```

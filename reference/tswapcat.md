# Trial-swap categorical randomization (tswapcat)

The trial-swap ("tswap") algorithm is a fixed–fixed randomization that
repeatedly attempts random 2×2 swaps until a valid one is found in each
iteration, reducing the number of wasted draws compared to the simple
swap. `tswapcat()` extends this logic to categorical matrices.

## Usage

``` r
tswapcat(
  x,
  n_iter = 1000L,
  output = c("category", "index"),
  swaps = "auto",
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

- seed:

  Integer used to seed random number generator, for reproducibility.

## Value

A matrix of the same dimensions as `x`, either randomized categorical
values (when `output = "category"`) or an integer index matrix
describing the permutation of entries (when `output = "index"`).

## References

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. *Ecology*, 81(9), 2606–2621.

Miklós, I. & Podani, J. (2004). Randomization of presence–absence
matrices: comments and new algorithms. *Ecology*, 85(1), 86–92.

Gotelli, N. J. & Entsminger, G. L. (2003). *EcoSim: Null models software
for ecology* (Version 7.0). Acquired Intelligence Inc. & Kesey-Bear,
Jericho (VT).

## See also

[`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md)
for an algorithm that produces equivalent results with better
computational efficiency.

## Examples

``` r
set.seed(123)
x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

# Randomize using swap algorithm
x_rand <- tswapcat(x, n_iter = 1000)

# Verify fixed-fixed constraint (row and column margins preserved)
all.equal(sort(x[1, ]), sort(x_rand[1, ]))
#> [1] TRUE
all.equal(sort(x[, 1]), sort(x_rand[, 1]))
#> [1] TRUE
```

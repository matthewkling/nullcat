# Categorical curveball randomization (curvecat)

Categorical generalization of the binary curveball algorithm (Strona et
al.) to matrices of categorical data. This function is a convenience
wrapper around
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
with `method = "curvecat"`.

## Usage

``` r
curvecat(x, n_iter = 1000L, output = "category", swaps = "auto", seed = NULL)
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

## Details

The curvecat algorithm randomizes a categorical matrix while keeping the
category multisets of each row and column fixed. In other words, the
permuted matrix has the same set of integer values in every row and
every column as the original matrix, but they are permuted. It operates
on pairs of rows at a time, grouping differing entries by unordered
category pairs and redistributing the orientation of those pairs while
preservingn the multiset of categories within each row. When there are
only two categories, `curvecat()` reduces to the behavior of the
original binary curveball algorithm applied to a 0/1 matrix.

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

# Use with a seed for reproducibility
x_rand1 <- curvecat(x, n_iter = 1000, seed = 42)
x_rand2 <- curvecat(x, n_iter = 1000, seed = 42)
identical(x_rand1, x_rand2)
#> [1] TRUE
```

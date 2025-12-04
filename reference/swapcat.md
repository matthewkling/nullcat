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

See also Gotelli & Entsminger (2003) *EcoSim: Null models software for
ecology* (Version 7.0) for implementation details of the binary swap
algorithm.

## See also

[`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md)
for an algorithm that produces equivalent results with better
computational efficiency.

## Examples

``` r
set.seed(123)
x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

# Randomize using swap algorithm
x_rand <- swapcat(x, n_iter = 1000)

# Verify fixed-fixed constraint (row and column margins preserved)
all.equal(sort(x[1, ]), sort(x_rand[1, ]))
#> [1] TRUE
all.equal(sort(x[, 1]), sort(x_rand[, 1]))
#> [1] TRUE
```

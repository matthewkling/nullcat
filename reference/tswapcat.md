# Trial-swap categorical randomization (tswapcat)

The trial-swap ("tswap") algorithm is a fixed–fixed randomization that
repeatedly attempts random 2×2 swaps until a valid one is found in each
iteration, reducing the number of wasted draws compared to the simple
swap. `tswapcat()` extends this logic to categorical matrices.

## Usage

``` r
tswapcat(x, n_iter = 1000L, output = c("category", "index"), seed = NULL)
```

## Arguments

- x:

  A matrix of categorical data, encoded as integers. Values should
  represent category or stratum membership for each cell.

- n_iter:

  Number of iterations. Default is 1000. Larger values yield more
  thorough mixing. Ignored for non-sequential methods. Minimum burn-in
  times can be estimated with
  [suggest_n_iter](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md).

- output:

  Character indicating type of result to return:

  - `"category"` (default) returns randomized matrix

  - `"index"` returns an index matrix describing where original entries
    moved.

- seed:

  Integer used to seed random number generator, for reproducibility.

## Value

A matrix of the same dimensions as `x`, either randomized categorical
values (when `output = "category"`) or an integer index matrix
describing the permutation of entries (when `output = "index"`).

## References

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. \*Ecology\*, 81(9), 2606–2621.

Miklós, I. & Podani, J. (2004). Randomization of presence–absence
matrices: comments and new algorithms. \*Ecology\*, 85(1), 86–92.

Gotelli, N. J. & Entsminger, G. L. (2003). \*EcoSim: Null models
software for ecology\* (Version 7.0). Acquired Intelligence Inc. &
Kesey-Bear, Jericho (VT).

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

# Column-constrained categorical randomization (c0cat)

\`r0cat()\` preserves the multiset of categories within each column but
randomizes their positions across rows, leaving row margins free. This
is the categorical analog to vegan's \`c0\` algorithm. It is a
non-sequential method.

## Usage

``` r
c0cat(x, n_iter = 1L, output = c("category", "index"), seed = NULL)
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

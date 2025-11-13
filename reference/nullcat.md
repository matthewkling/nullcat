# Categorical matrix randomization

Categorical generalizations of binary community null model algorithms.

## Usage

``` r
nullcat(
  x,
  method = nullcat_methods(),
  n_iter = 1000L,
  output = c("category", "index"),
  seed = NULL
)
```

## Arguments

- x:

  A matrix of categorical data, encoded as integers. Values should
  represent category or stratum membership for each cell.

- method:

  Character specifying the randomization algorithm to use. Options
  include:

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
    moved.

- seed:

  Integer used to seed random number generator, for reproducibility.

## Value

A matrix of the same dimensions as `x`, either randomized categorical
values (when `output = "category"`) or an integer index matrix
describing the permutation of entries (when `output = "index"`).

## Details

Note that categorical null models are for cell-level categorical data.
Site-level attributes (e.g., land cover) or species-level attributes
(e.g., functional traits) should be analyzed using different approaches.

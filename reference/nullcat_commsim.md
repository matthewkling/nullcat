# Nullcat-based commsim (non-sequential)

Construct a
[`vegan::commsim()`](https://vegandevs.github.io/vegan/reference/commsim.html)
object that uses
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
as a non-sequential null model for categorical / integer matrices. Each
simulated matrix is generated independently by applying
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
with `n_iter` trades starting from the original matrix.

## Usage

``` r
nullcat_commsim(
  n_iter = 10000,
  method = nullcat_methods(),
  output = c("category", "index"),
  row_groups = NULL,
  col_groups = NULL
)
```

## Arguments

- n_iter:

  Integer, number of iterations (trades) per simulated matrix. Must be a
  positive integer. Default is `1e4`.

- method:

  Character specifying which nullcat randomization algorithm to use. See
  [`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
  and
  [`nullcat_methods()`](https://matthewkling.github.io/nullcat/reference/nullcat_methods.md)
  for details.

- output:

  Character, passed to `nullcat(output = ...)`. Typically `"category"`
  (default) or `"index"`.

- row_groups, col_groups:

  Optional vectors of row or column groups, passed to
  [`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
  so that randomization happens independently within blocks. The data
  passed to
  [`vegan::nullmodel()`](https://vegandevs.github.io/vegan/reference/nullmodel.html)
  must have the matching number of rows and columns.

## Value

An object of class `"commsim"` suitable for use with
[`vegan::nullmodel()`](https://vegandevs.github.io/vegan/reference/nullmodel.html)
and
[`vegan::oecosimu()`](https://vegandevs.github.io/vegan/reference/oecosimu.html).

## Details

This generates a commsim object that is **non-sequential**: each
simulated matrix starts from the original matrix and is randomized
independently using `n_iter` trades of the chosen `method`.

When used via
[`vegan::simulate.nullmodel()`](https://vegandevs.github.io/vegan/reference/nullmodel.html),
the arguments behave as:

- `nsim`: number of simulated matrices to generate.

- `n_iter` (here, in `nullcat_commsim()`): number of trades per
  simulated matrix (controls how strongly each replicate is shuffled).

- `burnin` and `thin`: are **ignored** for this commsim, because
  `isSeq = FALSE` (the simulations are not a Markov chain).

In other words, treat `n_iter` as the tuning parameter for how
thoroughly each independent null matrix is randomized.

## See also

[`nullcat_batch()`](https://matthewkling.github.io/nullcat/reference/nullcat_batch.md)
if you just want a batch of null matrices without going through
**vegan**.

## Examples

``` r
# \donttest{
  library(vegan)
#> Loading required package: permute

  x  <- matrix(sample(1:5, 50, replace = TRUE), 10, 5)
  cs <- nullcat_commsim(n_iter = 1e4, method = "curvecat")

  nm   <- nullmodel(x, cs)
  sims <- simulate(nm, nsim = 999)
# }
```

# Categorical nullcat commsim (non-sequential)

Construct a
[`vegan::commsim`](https://vegandevs.github.io/vegan/reference/commsim.html)
object that uses
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
as a non-sequential null model for categorical / integer matrices. Each
simulated matrix is generated independently by applying
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
with `n_iter` trades starting from the original matrix.

## Usage

``` r
commsim_cat(
  n_iter = 10000,
  method = nullcat_methods(),
  output = c("category", "index")
)
```

## Arguments

- n_iter:

  Integer, number of iterations per simulated matrix.

- method:

  Character specifying which nullcat randomization algorithm to use. See
  [`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
  for details.

- output:

  Character, passed to `nullcat(output = ...)`. Typically `"category"`
  (default) or `"index"`.

## Value

An object of class `"commsim"` suitable for
[`vegan::nullmodel()`](https://vegandevs.github.io/vegan/reference/nullmodel.html).

# Bin quantitative data into strata

Bin quantitative data into strata

## Usage

``` r
stratify(
  x,
  breaks = NULL,
  n_strata = 5,
  transform = identity,
  offset = 0,
  zero_stratum = FALSE
)
```

## Arguments

- x:

  A matrix or vector containing non-negative values.

- breaks:

  Numeric vector of stratum breakpoints.

- n_strata:

  Integer giving the number of strata to split the data into. Must be 2
  or greater. Larger values yield randomizations with less mixing but
  higher fidelity to the original marginal distributions. Default is
  `5`. Ignored unless `breaks = NULL`.

- transform:

  A function used to transform the values in `x` before assigning them
  to `n_strata` equal-width intervals. Examples include `sqrt`, `log`,
  `rank` (for equal-occupancy strata), etc.; the default is `identity`.
  If `zero_stratum = TRUE`, the transformation is only applied on
  nonzero values. The function should pass NA values. This argument is
  ignored unless `breaks = NULL`.

- offset:

  Numeric value between -1 and 1 (default 0) indicating how much to
  shift stratum breakpoints relative to the binwidth (applied during
  quantization as: `breaks <- breaks + offset * bw`). To assess
  sensitivity to stratum boundaries, run
  [`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md)
  multiple times with different offset values. Ignored unless
  `breaks = NULL`.

- zero_stratum:

  Logical indicating whether to segregate zeros into their own stratum.
  If `FALSE` (the default), zeros will likely be combined into a stratum
  that also includes small positive numbers. If `breaks` is specified,
  zero simply gets added as an additional break; if not, one of the
  `n_strata` will represent zeros and the others will be nonzero ranges.

## Value

An object the same size as x, with integer values representing stratum
classifications.

## Examples

``` r
# Stratify a numeric vector
x <- c(0, 0, 0.1, 0.5, 1.2, 3.4, 5.6, 10.2)
stratify(x, n_strata = 3)
#> [1] 1 1 1 1 1 1 2 3
#> attr(,"breaks")
#> [1] -Inf  3.4  6.8  Inf

# With transformation
stratify(x, n_strata = 3, transform = log1p)
#> [1] 1 1 1 1 1 2 3 3
#> attr(,"breaks")
#> [1] -Inf  1.2  3.4  Inf

# Separate zero stratum
stratify(x, n_strata = 3, zero_stratum = TRUE)
#> [1] 1 1 2 2 2 2 3 3
#> attr(,"breaks")
#> [1] -Inf 0.00 5.15  Inf
```

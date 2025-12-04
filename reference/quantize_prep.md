# Prepare stratified null model overhead for quantize()

`quantize_prep()` precomputes all of the stratification and bookkeeping
needed by
[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md)
for a given quantitative community matrix. This is useful when you want
to generate many randomizations of the same dataset: the expensive steps
(strata assignment, value pools, and arguments for the underlying null
model) are computed once, and the resulting object can be passed to
`quantize(prep = ...)` for fast repeated draws.

## Usage

``` r
quantize_prep(
  x,
  method = nullcat_methods(),
  fixed = c("cell", "stratum", "row", "col"),
  breaks = NULL,
  n_strata = 5,
  transform = identity,
  offset = 0,
  zero_stratum = FALSE,
  n_iter = 1000
)
```

## Arguments

- x:

  Community matrix with sites in rows, species in columns, and
  nonnegative quantitative values in cells. This is the dataset for
  which stratification and null model overhead should be prepared.

- method:

  Character string specifying the null model algorithm. The default
  `"curvecat"` uses the categorical curveball algorithm. See
  [`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
  for alternative options.

- fixed:

  Character string specifying the level at which quantitative values are
  held fixed during randomization. One of:

  - `"cell"` (the default; only available when `method = "curvecat"`):
    values remain attached to their original cells and move with them
    during the categorical randomization. Row and column value
    distributions are not preserved, but the mapping between each
    original cell and its randomized destination is fixed.

  - `"stratum"`: values are shuffled globally within each stratum,
    holding only the overall stratum-level value distribution fixed.

  - `"row"`: values are shuffled within strata separately for each row,
    holding each row’s value multiset fixed. Not compatible with all
    `method`s.

  - `"col"`: values are shuffled within strata separately for each
    column, holding each column’s value multiset fixed.

  Note that this interacts with `method`: different null models fix
  different margins in the underlying binary representation.

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

- n_iter:

  Number of iterations. Default is 1000. Larger values yield more
  thorough mixing. Ignored for non-sequential methods. Minimum burn-in
  times can be estimated with
  [`suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md).

## Value

A list with class `"quantize_prep"` (if you want to set it) containing
the components needed by
[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md):

- `x`: original quantitative marix `x`,

- `strata`: integer matrix of the same dimension as `x`, giving the
  stratum index (`1:n_strata`) for each cell.

- `pool`: data structure encoding the quantitative value pools used
  during reassignment.

- `method`: the null model method used (as in the `method` argument).

- `n_strata`, `transform`, `offset`, `fixed`: the stratification and
  reassignment settings used to construct `strata` and `pool`.

- `sim_args`: named list of arguments passed to
  [`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
  (e.g. `n_iter`).

This object is intended to be passed unchanged to the `prep` argument of
[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md)
or
[`quantize_batch()`](https://matthewkling.github.io/nullcat/reference/quantize_batch.md).

## Details

Internally, `quantize_prep()`:

- transforms and stratifies `x` into `n_strata` numeric intervals (via
  [`stratify()`](https://matthewkling.github.io/nullcat/reference/stratify.md)),

- constructs the appropriate value pools given `fixed`, and

- assembles arguments for the underlying null model call to
  [`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md).

The returned object can be reused across calls to
[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md),
[`quantize_batch()`](https://matthewkling.github.io/nullcat/reference/quantize_batch.md),
or other helpers that accept a `prep` argument.

## See also

[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md),
[`quantize_batch()`](https://matthewkling.github.io/nullcat/reference/quantize_batch.md)

## Examples

``` r
set.seed(1)
comm <- matrix(rexp(50 * 40), nrow = 50,
               dimnames = list(paste0("site", 1:50),
                               paste0("sp", 1:40)))

# prepare overhead for a curvecat-backed stratified null model
prep <- quantize_prep(comm, method = "curvecat",
                      n_strata = 5,
                      fixed = "row",
                      n_iter = 2000)

# fast repeated randomizations using the same prep
rand1 <- quantize(prep = prep)
rand2 <- quantize(prep = prep)
```

# Stratified quantitative null models via quantize

`quantize()` implements a stratified randomization framework for
continuous ecological data. It discretizes quantitative values into
strata, randomizes the strata assignments using a categorical null model
algorithm (via
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)),
and then reassigns the original quantitative values according to the new
stratum layout.

## Usage

``` r
quantize(
  x = NULL,
  prep = NULL,
  method = nullcat_methods(),
  fixed = c("cell", "stratum", "row", "col"),
  breaks = NULL,
  n_strata = 5,
  transform = identity,
  offset = 0,
  zero_stratum = FALSE,
  n_iter = 1000,
  wt_row = NULL,
  wt_col = NULL,
  seed = NULL
)
```

## Arguments

- x:

  Community matrix with sites in rows, species in columns, and
  nonnegative quantitative values in cells. Can be `NULL` when `prep` is
  provided.

- prep:

  A `"quantize_prep"` object (from
  [`quantize_prep()`](https://matthewkling.github.io/nullcat/reference/quantize_prep.md)).
  If provided, `x` and all stratification arguments are ignored and the
  precomputed overhead is used directly for fast repeated draws.

- method:

  Character specifying the randomization algorithm to use. Options
  include the following; see details and linked functions for more info.

  - `"curvecat"`: categorical analog to `curveball`; see
    [`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md)
    for details.

  - `"swapcat"`: categorical analog to `swap`; see
    [`swapcat()`](https://matthewkling.github.io/nullcat/reference/swapcat.md)
    for details.

  - `"tswapcat"`: categorical analog to `tswap`; see
    [`tswapcat()`](https://matthewkling.github.io/nullcat/reference/tswapcat.md)
    for details.

  - `"r0cat"`: categorical analog to `r0`; see
    [`r0cat()`](https://matthewkling.github.io/nullcat/reference/r0cat.md)
    for details.

  - `"c0cat"`: categorical analog to `c0`; see
    [`c0cat()`](https://matthewkling.github.io/nullcat/reference/c0cat.md)
    for details.

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
    holding each row's value multiset fixed. Not compatible with all
    `method`s.

  - `"col"`: values are shuffled within strata separately for each
    column, holding each column's value multiset fixed.

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
  sensitivity to stratum boundaries, run `quantize()` multiple times
  with different offset values. Ignored unless `breaks = NULL`.

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

- wt_row:

  An optional square numeric matrix of non-negative weights controlling
  which pairs of rows are likely to exchange tokens during
  randomization. Must be `nrow(x)` by `nrow(x)`. This enables spatially
  or trait-constrained null models where nearby or similar sites
  exchange tokens more frequently.

  Values are treated as relative weights (not probabilities) and are
  normalized internally. The diagonal is ignored. The matrix should be
  symmetric. Only supported for sequential methods (`curvecat`,
  `swapcat`, `tswapcat`).

  When both `wt_row` and `wt_col` are supplied, `swaps` is forced to
  `"alternating"`, producing a Gibbs-like sweep that applies each weight
  matrix on its respective margin in alternation.

- wt_col:

  An optional square numeric matrix of non-negative weights controlling
  which pairs of columns are likely to exchange tokens during
  randomization. Must be `ncol(x)` by `ncol(x)`. See `wt_row` for
  details on weight interpretation.

- seed:

  Integer used to seed random number generator, for reproducibility.

## Value

A randomized version of `x`, with the same dimensions and dimnames. For
`method = "curvecat"`, the quantitative values are reassigned within
strata while preserving row and column stratum multisets. For binary
methods, the result corresponds to applying the chosen binary null model
to each stratum and recombining.

## See also

[`quantize_batch()`](https://matthewkling.github.io/nullcat/reference/quantize_batch.md)
for efficient generation of multiple randomized matrices;
[`quantize_commsim()`](https://matthewkling.github.io/nullcat/reference/quantize_commsim.md)
for integration with `vegan`.

## Examples

``` r
# toy quantitative community matrix
set.seed(1)
comm <- matrix(rexp(50 * 40), nrow = 50,
               dimnames = list(paste0("site", 1:50),
                               paste0("sp", 1:40)))

# default: curvecat-backed stratified randomization
rand1 <- quantize(comm)

# change stratification and preservation mode
rand2 <- quantize(comm, n_strata = 4,
                  transform = sqrt,
                  fixed  = "row",
                  n_iter    = 2000)

# use a different randomization algorithm
rand3 <- quantize(comm, method = "swapcat", n_iter = 10000)

# precompute overhead and reuse for many randomizations
prep  <- quantize_prep(comm, method = "curvecat",
                       n_strata = 5, fixed = "row")
rand4 <- quantize(prep = prep)
rand5 <- quantize(prep = prep)
```

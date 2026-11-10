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
  row_groups = NULL,
  col_groups = NULL,
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
    original cell and its randomized destination is fixed. Because cells
    only move when they trade places with a cell of a different stratum,
    a row or column whose values all fall in a single stratum is never
    changed;
    [`quantize_prep()`](https://matthewkling.github.io/nullcat/reference/quantize_prep.md)
    warns when this happens (see `zero_stratum` and `transform` for ways
    to separate low values into their own strata).

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
  When `row_groups` or `col_groups` is used, `n_iter` iterations are run
  separately on every block. The number needed can be more or fewer than
  for the full matrix, so estimate it by passing the same groups to
  [`suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md).

- row_groups, col_groups:

  Optional vectors assigning rows or columns of `x` to groups, so that
  randomization happens independently within each block (combination of
  a row group and a column group); see
  [`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md).
  Strata are defined once for the whole matrix, so every block uses the
  same stratum breaks, while the quantitative values are reassigned
  within each block (in every `fixed` mode, values never leave their
  block). With row groups representing regions, for example, each region
  keeps its own pool of species and values. `n_iter` iterations are run
  on every block.

- wt_row, wt_col:

  Deprecated. Experimental pairwise weight matrices previously used to
  bias which pairs of rows or columns exchange values. But weighting
  pair selection changes only how quickly a randomization mixes, not the
  null distribution it converges to; the only lasting effect of a weight
  matrix is to prevent exchanges between rows (or columns) that are not
  connected, directly or through intermediaries, by nonzero weights.
  Weight matrices are therefore now converted to `row_groups` /
  `col_groups` given by the connected components of their nonzero
  entries, with a warning. Use `row_groups` / `col_groups` instead.

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

# randomize within two regions of 25 sites each
region <- rep(1:2, each = 25)
rand6 <- quantize(comm, row_groups = region, n_iter = 2000)
all.equal(sort(comm[region == 1, ]), sort(rand6[region == 1, ]))
#> [1] TRUE
```

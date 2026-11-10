# Introduction to nullcat

``` r

library(nullcat)
```

## Introduction

`nullcat` provides null model algorithms for categorical and
quantitative community ecology data. It extends classic binary null
models (e.g., curveball, swap) to work with categorical data, and
introduces a stratified randomization framework for continuous data.

Categorical null models generalize binary null models to matrices where
cells contain integers representing discrete categories. Whereas binary
models preserve totals (e.g. row or column sums), categorical models
preserve **multisets** (e.g. category frequencies in each row or
column). For example, if a row contains \[1, 1, 2, 3, 3, 3\], a
row-preserving categorical algorithm would preserve this exact
multiset—two 1s, one 2, and three 3s—while randomizing where these
categories appear within the row.

Traditional binary models are simply the specific two-category cases of
their corresponding general categorical models (e.g. `curveball` is
equivalent to `curvecat` with only two categories).

### Use cases

Categorical null models are designed for **cell-level categorical
data**—characteristics of a species at a site—rather than site-level or
species-level attributes. They can be useful for categorical and ordinal
ecological data. Examples include:

- Ordinal states: phenological stages, health status, life stages
- Multi-state interactions: interaction type, infection status,
  mutualism outcomes
- Population genetic states: genotypes or allelic frequencies at the
  population level
- Binned quantitative measures: abundance classes, cover scales, or
  occurrence probability bins (the
  [`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md)
  function covered below automates binning and randomization)
- Binary data: species presence-absence data or similar

## Categorical null models

### Available algorithms

The package provides five categorical null model algorithms:

- **[`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md)**:
  Categorical curveball (most efficient)
- **[`swapcat()`](https://matthewkling.github.io/nullcat/reference/swapcat.md)**:
  Categorical 2×2 swap
- **[`tswapcat()`](https://matthewkling.github.io/nullcat/reference/tswapcat.md)**:
  Categorical trial-swap
- **[`r0cat()`](https://matthewkling.github.io/nullcat/reference/r0cat.md)**:
  Row-constrained randomization
- **[`c0cat()`](https://matthewkling.github.io/nullcat/reference/c0cat.md)**:
  Column-constrained randomization

### Example: Basic categorical randomization

``` r

# Create a categorical matrix
set.seed(123)
x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

# Randomize using curvecat (preserves row & column category multisets)
x_rand <- curvecat(x, n_iter = 1000)

# Verify margins are preserved
all.equal(sort(x[1, ]), sort(x_rand[1, ]))  # Row preserved
#> [1] TRUE
all.equal(sort(x[, 1]), sort(x_rand[, 1]))  # Column preserved
#> [1] TRUE
```

### Understanding fixed margins

Different algorithms preserve different margins:

``` r

set.seed(456)
x <- matrix(sample(1:3, 60, replace = TRUE), nrow = 10)

# Fixed-fixed: both row and column margins preserved
x_fixed <- curvecat(x, n_iter = 1000)

# Row-constrained: only row margins preserved
x_row <- r0cat(x)

# Column-constrained: only column margins preserved
x_col <- c0cat(x)

# Check row sums by category
table(x[1, ])           # Original row 1
#> 
#> 1 2 3 
#> 2 2 2
table(x_fixed[1, ])     # Preserved in curvecat
#> 
#> 1 2 3 
#> 2 2 2
table(x_row[1, ])       # Preserved in r0cat
#> 
#> 1 2 3 
#> 2 2 2
table(x_col[1, ])       # Not preserved in c0cat
#> 
#> 1 2 3 
#> 1 3 2
```

## Quantitative null models

The
[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md)
function provides stratified randomization for continuous community data
(abundance, biomass, cover, occurrence probability). It works by
applying the following steps:

1.  **Stratify**: use a flexible binning scheme to convert continuous
    values into discrete strata
2.  **Randomize**: apply one of the categorical null models described
    above to the stratified matrix
3.  **Reassign**: map randomized strata back to original values, subject
    to a chosen constraint

**Why use quantize?** Existing quantitative null models lack flexibility
in what margins they preserve. While binary curveball elegantly
preserves row and column totals, and vegan’s quasiswap models do this
for count data, there’s no framework for preserving row and/or column
**value distributions** for continuous data (abundance, biomass, cover).
The stratification approach fills this gap: when using `fixed = "row"`
or `fixed = "col"`, one dimension’s value multisets are preserved
exactly while the other is preserved at the resolution of strata. This
is the closest available approximation to a quantitative fixed-fixed
null model. Additionally, stratification provides explicit control over
the preservation-randomization tradeoff through the number of strata
used.

### Example: Basic quantitative randomization

``` r

# Create a quantitative community matrix
set.seed(789)
comm <- matrix(runif(200), nrow = 20)

# Default: curvecat-backed stratified randomization with 5 strata
rand1 <- quantize(comm, n_strata = 5, n_iter = 2000, fixed = "row")

# Values are similar but rearranged
cor(as.vector(comm), as.vector(rand1))
#> [1] 0.1757922
plot(rowSums(comm), rowSums(rand1))
```

![](nullcat_files/figure-html/quantitative-basic-1.png)

``` r

plot(colSums(comm), colSums(rand1))
```

![](nullcat_files/figure-html/quantitative-basic-2.png)

### Stratification options

When you call
[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md),
it uses the helper function
[`stratify()`](https://matthewkling.github.io/nullcat/reference/stratify.md)
to convert your continuous data to discrete bins. The stratification
scheme can strongly influence your randomization results, because it is
frequencies of these categories that are constrained during matrix
permutation. By default, data are split into a handful of equal-width
bins, but this can be controlled with parameters including `n_strata`,
`transform`, `zero_stratum`, `breaks`, and `offset`. Let’s call
[`stratify()`](https://matthewkling.github.io/nullcat/reference/stratify.md)
directly here, to understand its behavior—a good idea for a real
analysis, not just a tutorial.

``` r

set.seed(200)
x <- rexp(100, .1)

# More strata = less mixing, higher fidelity to original distribution
s3 <- stratify(x, n_strata = 3)
s10 <- stratify(x, n_strata = 10)

table(s3)   # Coarser bins
#> s3
#>  1  2  3 
#> 89  9  2
table(s10)  # Finer bins
#> s10
#>  1  2  3  4  5  6 10 
#> 49 22 17  4  4  2  2

# Transform before stratifying (e.g., log-transform for skewed data)
s_log <- stratify(x, n_strata = 5, transform = log1p)
table(s_log)
#> s_log
#>  1  2  3  4  5 
#> 14 26 27 27  6

# Rank transform creates equal-occupancy strata
s_rank <- stratify(x, n_strata = 5, transform = rank)
table(s_rank)  # Nearly equal counts per stratum
#> s_rank
#>  1  2  3  4  5 
#> 20 20 20 20 20

# Separate zeros into their own stratum
x_with_zeros <- c(0, 0, 0, x)
s_zero <- stratify(x_with_zeros, n_strata = 4, zero_stratum = TRUE)
table(s_zero)
#> s_zero
#>  1  2  3  4 
#>  3 89  9  2
```

### Controlling what’s preserved

While the mechanics of the chosen randomization algorithm determine
which margins of the randomized *categories* are preserved, the `fixed`
parameter controls how the original *continuous* values are mapped back
onto the randomized categories, determining which margins preserve their
original value distributions. There are four options for the `fixed`
parameter: row, col, cell, and stratum, though only a subset of these
are relevant for some algorithms. For a method like the default
`curvecat` algorithm that preserves both row and column categorical
multisets, all four are available:

``` r

set.seed(100)
comm <- matrix(rexp(100), nrow = 10)

# Preserve row value multisets (quantitative row sums maintained)
rand_row <- quantize(comm, n_strata = 5, fixed = "row", n_iter = 2000)
all.equal(rowSums(comm), rowSums(rand_row))
#> [1] TRUE

# Preserve column value multisets (quantitative column sums maintained)
rand_col <- quantize(comm, n_strata = 5, fixed = "col", n_iter = 2000)
all.equal(colSums(comm), colSums(rand_col))
#> [1] TRUE

# Cell-level preservation: each value moves with its original cell location
# The categorical randomization determines WHERE cells go, but each cell
# carries its original value with it
rand_cell <- quantize(comm, n_strata = 5, fixed = "cell", n_iter = 2000)

# Values shuffled globally within strata, holding none of the above fixed
rand <- quantize(comm, n_strata = 5, fixed = "stratum", n_iter = 2000)

# For non-fixed-fixed methods like r0cat or c0cat, only some fixed options make sense:
# r0cat with fixed="col" or c0cat with fixed="row" would be incompatible
```

## Efficient repeated randomization

The quantize routine includes some setup, including stratification, that
doesn’t need to be repeatedly recomputed when generating multiple
randomized versions of the same matrix. For generating null
distributions, you can pre-compute overhead once using
[`quantize_prep()`](https://matthewkling.github.io/nullcat/reference/quantize_prep.md),
or you can use `*_batch()` functions to efficiently generate multiple
replicates (optionally using parallel processing).

``` r

set.seed(400)
comm <- matrix(rexp(200), nrow = 20)

# Prepare once
prep <- quantize_prep(comm, n_strata = 5, fixed = "row", n_iter = 2000)

# Generate many randomizations efficiently
rand1 <- quantize(prep = prep)
rand2 <- quantize(prep = prep)
rand3 <- quantize(prep = prep)

# Or use batch functions
nulls <- quantize_batch(comm, n_reps = 99, n_strata = 5, 
                       fixed = "row", n_iter = 2000)
dim(nulls)  # 20 rows × 10 cols × 99 replicates
#> [1] 20 10 99
```

## Convergence diagnostics

Sequential algorithms (curvecat, swapcat, tswapcat) require sufficient
iterations to reach stationarity. You can use the
[`trace_cat()`](https://matthewkling.github.io/nullcat/reference/trace_cat.md)
function with either `nullcat` or `quantize` to assess how many
iterations are needed for your particular dataset and algorithm. The
trace shows how similarity to the original matrix changes over
iterations; when the trace stabilizes around zero, the algorithm has
reached stationarity. You can also use
[`suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md)
to get an automated suggestion for a reasonable number of burn-in
iterations to reach stationarity.

``` r

set.seed(300)
x <- matrix(sample(1:5, 400, replace = TRUE), nrow = 20)

# Generate trace showing mixing over iterations
trace <- trace_cat(x, fun = "nullcat", method = "curvecat",
                   n_iter = 1000, n_chains = 3, thin = 10)

# Visual inspection
plot(trace)
```

![](nullcat_files/figure-html/trace-diagnostics-1.png)

``` r


# Automatic burn-in suggestion
suggested <- suggest_n_iter(trace, tail_frac = 0.3)
print(suggested)
#> suggested_n_iter object
#> -----------------------
#> Converged: TRUE 
#> Suggested n iterations: 410
```

## Constrained randomization within groups

Standard null models let any row exchange tokens with any other row. In
many ecological contexts this is unrealistic: for example, species
assemblages at sites 3000 km apart shouldn’t contribute equally to a
null expectation designed to test *local* community assembly. `nullcat`
handles this by randomizing within groups. The `row_groups` and
`col_groups` arguments assign rows (e.g. sites) and/or columns
(e.g. species) to groups; the matrix is split into blocks (every
combination of a row group and a column group), and each block is
randomized independently. Tokens therefore never leave their block.

For algorithms that fix row and column margins, this adds a constraint
to the null model. With sites grouped into regions, for example, each
region keeps the category counts of every species: a null model with a
regional species pool. Grouping works with every `nullcat` method, with
[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md),
and with the batch, trace, and `vegan` integration functions.

### Fixed groups

``` r

set.seed(42)

# 20 sites x 15 species categorical matrix
x <- matrix(sample(1:4, 300, replace = TRUE), nrow = 20)

# Sites belong to two regions
region <- rep(c("north", "south"), each = 10)

x_reg <- curvecat(x, n_iter = 1000, row_groups = region)

# Row and column margins are still preserved...
all.equal(sort(x[1, ]), sort(x_reg[1, ]))
#> [1] TRUE
all.equal(sort(x[, 1]), sort(x_reg[, 1]))
#> [1] TRUE

# ...and so are the category counts of each species within each region
all.equal(rowsum((x == 1) * 1, region), rowsum((x_reg == 1) * 1, region))
#> [1] TRUE
```

Column groups work the same way, e.g. restricting exchanges to species
within the same clade or guild. Supplying both gives independent
randomization within every region-by-clade block.

When groups are used, `n_iter` iterations are run on every block. That
can be more or fewer iterations than the whole matrix needs: a smaller
block has fewer rows to mix among, but also fewer possible swaps. Given
the same groups,
[`trace_cat()`](https://matthewkling.github.io/nullcat/reference/trace_cat.md)
and
[`suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md)
trace the grouped randomization as it will be run and suggest the number
of iterations per block:

``` r

suggest_n_iter(x = x, n_iter = 2000, row_groups = region)
```

### Random spatial tilings

Fixed regions impose arbitrary boundaries: two neighboring sites on
either side of a boundary never exchange tokens, while two sites at
opposite ends of a region always can. A smoother, distance-based
alternative is to use a different random grouping in each replicate.
[`tile_groups()`](https://matthewkling.github.io/nullcat/reference/tile_groups.md)
covers the sites with randomly placed and oriented tiles, returning one
tiling per column. Batch functions accept such a matrix of groups and
use column `r` for replicate `r`:

``` r

set.seed(42)

# 200 sites with coordinates, 15 species
coords <- cbind(runif(200, 0, 100), runif(200, 0, 100))
x <- matrix(sample(1:4, 3000, replace = TRUE), nrow = 200)

# one random tiling per replicate, with tiles 30 units apart
tiles <- tile_groups(coords, spacing = 30, n = 99, lonlat = FALSE)

sims <- nullcat_batch(x, n_reps = 99, row_groups = tiles, n_iter = 500)
dim(sims)
#> [1] 200  15  99
```

Across tilings, two sites a distance `d` apart share a tile with a
probability that falls from 1 at `d = 0` to about 0.5 at
`d = 0.42 * spacing`, and to 0 at about `1.15 * spacing`. This acts as
the distance kernel of the null model: each site is pooled with its
neighbors, never with distant sites. A vector of spacings (one per
tiling) blends tile sizes for a smoother, longer-tailed kernel.

The example above uses planar coordinates (`lonlat = FALSE`), which are
tiled with regular hexagons. For longitude/latitude coordinates
(`lonlat = TRUE`, the default), tiles are instead laid out on the
sphere, so results don’t depend on a map projection, and `spacing` is in
meters. Spherical tiles all have the same area as a hexagon with the
given spacing, but their shapes vary between irregular hexagons and
near-squares, depending on where sites fall on the underlying Fibonacci
lattice. Here are six random tilings of a regular grid of sites, with
tiles colored arbitrarily:

``` r

set.seed(123)
s <- seq(0, 10, by = 0.25)
ll <- expand.grid(lon = s, lat = s)
tiles_ll <- tile_groups(ll, spacing = 300e3, n = 400)

par(mfrow = c(2, 3), mar = c(0.5, 0.5, 0.5, 0.5))
for (i in 1:6) {
  image(s, s, matrix(tiles_ll[, i] %% 7, length(s)), col = hcl.colors(7, "viridis"),
        asp = 1, axes = FALSE, xlab = "", ylab = "")
}
```

![](nullcat_files/figure-html/groups-tiles-sphere-1.png)

Averaged over many tilings, the varying shapes and orientations wash
out, giving practically the same smooth, circular distance kernel as
regular hexagons (with a slightly longer tail, reaching 0 at about
`1.3 * spacing`). Here is the proportion of tilings in which each site
shares a tile with the site at the center of the grid:

``` r

focal <- which.min((ll$lon - 5)^2 + (ll$lat - 5)^2)
p_shared <- colMeans(t(tiles_ll) == tiles_ll[focal, ])

filled.contour(s, s, matrix(p_shared, length(s)), asp = 1,
               color.palette = function(n) hcl.colors(n, "viridis"),
               plot.title = title(xlab = "longitude", ylab = "latitude"),
               key.title = title(main = "shared\ntile", cex.main = 0.8))
```

![](nullcat_files/figure-html/groups-kernel-1.png)

Sites near edges or gaps tend to fall in small tiles, and a site alone
in its tile cannot be changed by fixed-margin algorithms.
[`group_sizes()`](https://matthewkling.github.io/nullcat/reference/group_sizes.md)
gives each site’s average pool size, which is worth checking before
choosing `spacing`:

``` r

summary(group_sizes(tiles))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    3.98   10.98   13.58   13.52   16.07   19.60
```

## Integration with `vegan`

`nullcat` includes several functions to integrate with the `vegan`
package’s null model framework. For both
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
and
[`quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.md),
there are sequential and non-sequential variants of these functions that
create
[`vegan::commsim()`](https://vegandevs.github.io/vegan/reference/commsim.html)
objects suitable for use with
[`vegan::nullmodel()`](https://vegandevs.github.io/vegan/reference/nullmodel.html)
and
[`vegan::oecosimu()`](https://vegandevs.github.io/vegan/reference/oecosimu.html).
Examples include:

``` r

library(vegan)

# Categorical data
x_cat <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
cs <- nullcat_commsim_seq(method = "curvecat")
nm <- nullmodel(x_cat, cs)
sims <- simulate(nm, nsim = 99, burnin = 1000, thin = 100)

# Quantitative data
x_quant <- matrix(rexp(100), nrow = 10)
cs_quant <- quantize_commsim(n_strata = 5, method = "curvecat", n_iter = 2000)
nm_quant <- nullmodel(x_quant, cs_quant)
sims_quant <- simulate(nm_quant, nsim = 99)
```

## References

- Strona et al. (2014). A fast and unbiased procedure to randomize
  ecological binary matrices. *Nature Communications*, 5, 4114.
- Gotelli (2000). Null model analysis of species co-occurrence patterns.
  *Ecology*, 81(9), 2606-2621.

## Further reading

- Package website: <https://matthewkling.github.io/nullcat/>
- Report issues: <https://github.com/matthewkling/nullcat/issues>

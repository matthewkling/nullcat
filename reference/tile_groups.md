# Random hexagonal tilings for spatially constrained null models

Generates random groupings of sites by tiling space with regular
hexagons that are randomly positioned and oriented, independently for
each column of the result. Passing the result as `row_groups` to
[`nullcat_batch()`](https://matthewkling.github.io/nullcat/reference/nullcat_batch.md)
or
[`quantize_batch()`](https://matthewkling.github.io/nullcat/reference/quantize_batch.md)
randomizes each replicate within the tiles of its own tiling, so that
sites only exchange values with nearby sites, and the resulting null
distribution blends over many tile placements.

## Usage

``` r
tile_groups(coords, spacing, n = 1L, lonlat = TRUE)
```

## Arguments

- coords:

  A two-column numeric matrix (or data frame) of site coordinates, one
  row per row of the community matrix. With `lonlat = TRUE`, longitude
  and latitude in decimal degrees; otherwise planar (projected) x and y
  coordinates.

- spacing:

  Distance between the centres of neighbouring tiles. With
  `lonlat = TRUE` this is in metres (great-circle distance on a
  spherical Earth of radius 6,371,008.8 m); otherwise it is in the units
  of `coords`. Either a single value, or a vector of length `n` giving
  the spacing for each tiling (e.g. random draws, for a blend of
  scales).

- n:

  Number of tilings (columns of the result). Use one per randomization
  replicate.

- lonlat:

  Logical: are `coords` longitude/latitude (`TRUE`, the default) or
  planar coordinates (`FALSE`)?

## Value

An integer matrix with one row per site and `n` columns. Each column is
one tiling, giving the tile membership of every site (tile labels are
arbitrary and only meaningful within a column).

## Details

Each tiling is built from a set of regularly spaced seed points, and
every site is assigned to its nearest seed. With `lonlat = TRUE`, seeds
form a Fibonacci lattice on the sphere with a uniformly random 3-D
rotation for each tiling. This gives near-hexagonal, near-equal-area
tiles of any size, anywhere on the globe, with no projection distortion;
nearest seeds are found by a fast exact search, so fine spacings are no
more costly than coarse ones. With `lonlat = FALSE`, seeds form a planar
hexagonal grid with a random rotation and offset for each tiling.

Because tile placement is random, two sites a distance `d` apart share a
tile in some tilings and not others. Across tilings, the probability
that they share a tile falls from 1 at `d = 0` to about 0.5 at
`d = 0.42 * spacing`, and to 0 at about the tile diameter:
`1.15 * spacing` for a regular hexagon (the planar grid). Tiles of the
spherical lattice all have nearly the same area, but some are moderately
elongated, so a few pairs of sites up to about `1.3 * spacing` apart can
occasionally share a tile. This "co-tile" probability acts as the
distance kernel of the null model: sites are pooled with their
neighbours, and never with distant sites. Supplying a vector of spacings
blends tiles of several sizes, which smooths and lengthens the kernel.

Sites near edges, gaps, or coastlines tend to fall in tiles with few
other sites. A site alone in its tile cannot be changed by fixed-margin
randomizations in that replicate. Use
[`group_sizes()`](https://matthewkling.github.io/nullcat/reference/group_sizes.md)
to check the average pool size of each site before running
randomizations. When randomizing within groups, `n_iter` is applied to
every tile, so it only needs to be large enough for the largest tile;
[`suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md)
given the same groups traces the largest tile.

## See also

[`group_sizes()`](https://matthewkling.github.io/nullcat/reference/group_sizes.md),
[`nullcat_batch()`](https://matthewkling.github.io/nullcat/reference/nullcat_batch.md),
[`quantize_batch()`](https://matthewkling.github.io/nullcat/reference/quantize_batch.md)

## Examples

``` r
# planar coordinates
set.seed(1)
xy <- cbind(runif(300, 0, 100), runif(300, 0, 100))
tiles <- tile_groups(xy, spacing = 30, n = 99, lonlat = FALSE)
summary(group_sizes(tiles))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   5.909  17.876  21.985  20.816  24.235  29.485 

# randomize a community matrix within a different tiling in each replicate
comm <- matrix(rbinom(300 * 20, 1, 0.3), 300)
sims <- nullcat_batch(comm, n_reps = 99, row_groups = tiles, n_iter = 500)

# geographic coordinates, with 200 km spacing
ll <- cbind(runif(300, -124, -114), runif(300, 32, 42))
tiles_ll <- tile_groups(ll, spacing = 200e3, n = 99)
```

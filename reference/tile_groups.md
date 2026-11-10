# Random spatial tilings for spatially constrained null models

Generates random groupings of sites by partitioning space into compact,
equal-area tiles that are randomly positioned and oriented,
independently for each column of the result. Passing the result as
`row_groups` to
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

  The size of the tiles. With `lonlat = FALSE`, the distance between the
  centers of neighboring hexagonal tiles, in the units of `coords`. With
  `lonlat = TRUE`, the equivalent spacing in meters: tiles have the same
  area as hexagons with this center spacing (`sqrt(3) / 2 * spacing^2`),
  on a spherical Earth of radius 6,371,008.8 m. Either a single value,
  or a vector of length `n` giving the spacing for each tiling (e.g.
  random draws, for a blend of scales).

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

Each tiling is built from a set of evenly spaced seed points, and every
site is assigned to its nearest seed, so tiles are the Voronoi cells of
the seeds.

With `lonlat = FALSE`, seeds form a regular hexagonal grid with a random
rotation and offset for each tiling, so tiles are regular hexagons.

With `lonlat = TRUE`, seeds form a Fibonacci lattice on the sphere with
a uniformly random 3-D rotation for each tiling. This works at any
extent, up to the whole globe, with no projection distortion, and
nearest seeds are found by a fast exact search, so fine spacings are no
more costly than coarse ones. Tiles all have the same area, but not the
same shape: the lattice is only approximately hexagonal, so depending on
where sites fall on it, tiles range from irregular hexagons to
near-squares.

Because tile placement is random, two sites a distance `d` apart share a
tile in some tilings and not others. Since every tiling is placed and
oriented at random, the probability that they share a tile depends only
on `d`, not on where the sites are or the direction between them. It
falls from 1 at `d = 0` to about 0.5 at `d = 0.42 * spacing` and about
0.17 at `d = 0.75 * spacing`, and is near zero beyond `d = spacing`. It
reaches zero at the tile diameter: `1.15 * spacing` for regular
hexagons, or about `1.3 * spacing` on the sphere, where some tiles are
slightly elongated. Otherwise the two give practically the same curve.
This "co-tile" probability acts as the distance kernel of the null
model: sites are pooled with their neighbors, and never with distant
sites. Supplying a vector of spacings blends tiles of several sizes,
which smooths and lengthens the kernel.

Sites near edges, gaps, or coastlines tend to fall in tiles with few
other sites. A site alone in its tile cannot be changed by fixed-margin
randomizations in that replicate. Use
[`group_sizes()`](https://matthewkling.github.io/nullcat/reference/group_sizes.md)
to check the average pool size of each site before running
randomizations. When randomizing within groups, `n_iter` is applied to
every tile; estimate it by passing the same groups to
[`suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md),
since tiles can need more or fewer iterations than the full matrix.

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

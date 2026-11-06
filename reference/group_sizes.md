# Mean group size for each row

For a grouping vector or a matrix of groupings (such as the output of
[`tile_groups()`](https://matthewkling.github.io/nullcat/reference/tile_groups.md),
with one column per randomization), returns the average number of
members in each element's group. With a groups matrix, this is the
average over columns. It describes the size of the pool each row (e.g.
site) is randomized within, and depends only on the grouping, so it can
be checked before running any randomizations, for example while choosing
a tile `spacing`. Rows that are often in very small groups are rarely
changed by a fixed-margin randomization; a row alone in its group is
never changed.

## Usage

``` r
group_sizes(groups)
```

## Arguments

- groups:

  A vector of group labels, or a matrix of group labels with one column
  per randomization.

## Value

A numeric vector with one value per element (row) of `groups`.

## See also

[`tile_groups()`](https://matthewkling.github.io/nullcat/reference/tile_groups.md)

## Examples

``` r
group_sizes(c("a", "a", "b", "c", "c", "c"))
#> [1] 2 2 1 3 3 3

set.seed(1)
coords <- cbind(runif(200, 0, 10), runif(200, 0, 10))
g <- tile_groups(coords, spacing = 3, n = 50, lonlat = FALSE)
summary(group_sizes(g))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    2.88   12.46   15.18   15.05   18.02   22.72 
```

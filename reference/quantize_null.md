# Generate a null distribution using quantize()

Generate a null distribution using quantize()

## Usage

``` r
quantize_null(x, n_reps = 999L, stat = NULL, n_cores = 1L, ...)
```

## Arguments

- x:

  Community matrix (species × sites, or any numeric matrix).

- n_reps:

  Number of randomizations to generate. Default is \`999\`.

- stat:

  Optional summary function taking a matrix and returning a numeric
  statistic (e.g. \`rowSums\` with abundance data would give total
  abundance per site). If \`NULL\` (default), the function returns the
  full set of randomized matrices.

- n_cores:

  Number of compute cores to use for parallel processing. Default is
  \`1\`.

- ...:

  Additional arguments passed to \`quantize()\` (e.g. \`method\`,
  \`breaks\`, \`n_strata\`, \`transform\`, \`offset\`, \`zero_stratum\`,
  \`fixed\`, \`n_iter\`, etc.).

## Value

If stat is NULL: a 3D array (rows × cols × n_reps). If stat is not NULL:
a numeric array of statistic values (dimensionality will depend on
stat).

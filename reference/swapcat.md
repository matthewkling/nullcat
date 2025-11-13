# Categorical swap randomization (swapcat)

Categorical generalization of the binary 2x2 swap algorithm to matrices
of categorical data. This function is a convenience wrapper around
\[nullcat()\] with \`method = "swapcat"\`.

## Usage

``` r
swapcat(x, n_iter = 1000L, output = c("category", "index"), seed = NULL)
```

## Arguments

- x:

  A matrix of categorical data, encoded as integers. Values should
  represent category or stratum membership for each cell.

- n_iter:

  Number of iterations. Default is 1000. Larger values yield more
  thorough mixing. Ignored for non-sequential methods. Minimum burn-in
  times can be estimated with
  [suggest_n_iter](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.md).

- output:

  Character indicating type of result to return:

  - `"category"` (default) returns randomized matrix

  - `"index"` returns an index matrix describing where original entries
    moved.

- seed:

  Integer used to seed random number generator, for reproducibility.

## Value

A matrix of the same dimensions as `x`, either randomized categorical
values (when `output = "category"`) or an integer index matrix
describing the permutation of entries (when `output = "index"`).

## Details

The swapcat algorithm attempts random 2x2 swaps of the form

\$\$ \begin{pmatrix} a & b \\ b & a \end{pmatrix} \leftrightarrow
\begin{pmatrix} b & a \\ a & b \end{pmatrix} \$\$

where \\a\\ and \\b\\ are distinct categories. These swaps preserve the
multiset of categories in each row and column. With only two categories
present, \`swapcat()\` reduces to the behavior of the standard binary
swap algorithm.

## References

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. \*Ecology\*, 81(9), 2606–2621.

See also Gotelli & Entsminger (2003) \*EcoSim: Null models software for
ecology\* (Version 7.0) for implementation details of the binary swap
algorithm.

## See also

\[nullcat()\], \[nullcat_methods()\]

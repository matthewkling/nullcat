# Nullcat-based commsim (non-sequential)

Construct a \`vegan::commsim()\` object that uses \[nullcat()\] as a
non-sequential null model for categorical / integer matrices. Each
simulated matrix is generated independently by applying \[nullcat()\]
with \`n_iter\` trades starting from the original matrix.

## Usage

``` r
nullcat_commsim(
  n_iter = 10000,
  method = nullcat_methods(),
  output = c("category", "index")
)
```

## Arguments

- n_iter:

  Integer, number of iterations (trades) per simulated matrix. Must be a
  positive integer. Default is \`1e4\`.

- method:

  Character specifying which nullcat randomization algorithm to use. See
  \[nullcat()\] and \[nullcat_methods()\] for details.

- output:

  Character, passed to \`nullcat(output = ...)\`. Typically
  \`"category"\` (default) or \`"index"\`.

## Value

An object of class \`"commsim"\` suitable for use with
\`vegan::nullmodel()\` and \`vegan::oecosimu()\`.

## Details

This generates a commsim object that is \*\*non-sequential\*\*: each
simulated matrix starts from the original matrix and is randomized
independently using \`n_iter\` trades of the chosen \`method\`.

When used via \`vegan::simulate.nullmodel()\`, the arguments behave as:

- \`nsim\`: number of simulated matrices to generate.

- \`n_iter\` (here, in \`nullcat_commsim()\`): number of trades per
  simulated matrix (controls how strongly each replicate is shuffled).

- \`burnin\` and \`thin\`: are \*\*ignored\*\* for this commsim, because
  \`isSeq = FALSE\` (the simulations are not a Markov chain).

In other words, treat \`n_iter\` as the tuning parameter for how
thoroughly each independent null matrix is randomized.

## See also

\[nullcat_batch()\] if you just want a batch of null matrices without
going through \*\*vegan\*\*.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(vegan)

  x  <- matrix(sample(1:5, 50, replace = TRUE), 10, 5)
  cs <- nullcat_commsim(n_iter = 1e4, method = "curvecat")

  nm   <- nullmodel(x, cs)
  sims <- simulate(nm, nsim = 999)
} # }
```

# Quantize-based commsim (non-sequential)

Construct a \[vegan::commsim()\] object that uses \[quantize()\] as a
non-sequential null model for numeric community matrices. Each simulated
matrix is generated independently by applying \[quantize()\] with
\`n_iter\` trades (via its internal call to \[nullcat()\]) starting from
the original matrix.

## Usage

``` r
quantize_commsim(n_iter = 10000, ...)
```

## Arguments

- n_iter:

  Integer, number of iterations (trades) per simulated matrix. Must be a
  positive integer. Default is \`1e4\`.

- ...:

  Arguments passed to \[quantize()\], such as \`breaks\`, \`n_strata\`,
  \`transform\`, \`offset\`, \`zero_stratum\`, \`fixed\`, \`method\`,
  etc. Do \*\*not\*\* supply \`x\` or \`n_iter\` here; these are set
  internally by \`quantize_commsim()\`. See \[quantize()\] for details.

## Value

An object of class \`"commsim"\` suitable for \[vegan::nullmodel()\] and
\[vegan::oecosimu()\].

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

\[quantize_batch()\] if you just want a batch of null matrices without
going through \*\*vegan\*\*.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(vegan)

  x <- matrix(rexp(50), 10, 5)

  cs <- quantize_commsim(
    n_strata = 10,
    method   = "curvecat",
    n_iter   = 1000L
  )

  nm   <- nullmodel(x, cs)
  sims <- simulate(nm, nsim = 999)
} # }
```

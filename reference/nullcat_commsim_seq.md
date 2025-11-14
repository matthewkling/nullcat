# Nullcat-based commsim (sequential / Markov chain)

Construct a \`vegan::commsim()\` object that uses \[nullcat()\] as a
\*sequential\* null model: successive simulated matrices form a Markov
chain. Internally, each simulation "step" advances the chain by \`thin\`
trades of the chosen \`method\` (e.g. \`"curvecat"\`), where \`thin\` is
supplied via \`vegan::simulate.nullmodel()“ arguments. This is analogous
to how sequential swap / curveball null models are used in
\*\*vegan\*\*, but extended to categorical data via \[nullcat()\].

## Usage

``` r
nullcat_commsim_seq(
  method = nullcat_methods(),
  output = c("category", "index")
)
```

## Arguments

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

This model is \*\*sequential\*\*: simulated matrices form a Markov
chain. The current matrix is updated in-place by repeated calls to the
randomization model, and successive matrices are obtained by advancing
the chain.

In \`vegan::simulate.nullmodel()\`, the control arguments behave as:

- \`nsim\`: number of matrices to \*store\* from the chain.

- \`thin\`: number of trades per step. Each "step" of the chain applies
  \`thin\` trades of the chosen \`method\` to the current state before
  possibly storing it.

- \`burnin\`: number of initial steps to perform (each with \`thin\`
  trades) before storing any matrices, i.e. the Markov chain burn-in.

There is no \`n_iter\` argument here: mixing is controlled entirely by
\`thin\` (trades per step) and \`burnin\` (number of initial steps
discarded), in the same spirit as sequential swap / curveball models in
\*\*vegan\*\*.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(vegan)

  x  <- matrix(sample(1:5, 50, replace = TRUE), 10, 5)
  cs <- nullcat_commsim_seq(method = "curvecat")

  nm <- nullmodel(x, cs)

  # control the chain with 'thin' and 'burnin'
  sims <- simulate(nm, nsim = 999, thin = 100, burnin = 1000)
} # }
```

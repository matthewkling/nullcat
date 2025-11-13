# Categorical curveball commsim (sequential / Markov chain)

Construct a
[`vegan::commsim`](https://vegandevs.github.io/vegan/reference/commsim.html)
object that uses
[`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
as a \*sequential\* null model: successive simulated matrices form a
Markov chain, and mixing is controlled via the `thin` and `burnin`
arguments to `simulate.nullmodel()`.

## Usage

``` r
commsim_cat_seq(method = nullcat_methods(), output = c("category", "index"))
```

## Arguments

- method:

  Character specifying which nullcat randomization algorithm to use. See
  [`nullcat()`](https://matthewkling.github.io/nullcat/reference/nullcat.md)
  for details.

- output:

  Character, passed to `nullcat(output = ...)`. Typically `"category"`
  (default) or `"index"`.

## Value

An object of class `"commsim"` suitable for
[`vegan::nullmodel()`](https://vegandevs.github.io/vegan/reference/nullmodel.html).

## Details

Internally, each simulation "step" advances the chain by `thin` curvecat
trades.

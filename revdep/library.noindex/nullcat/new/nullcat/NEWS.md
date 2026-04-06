# nullcat 0.2.0

* New feature: weighted pair sampling via `wt_row` and `wt_col` parameters.
  Enables spatially, phylogenetically, or trait-constrained null models by
  weighting which pairs of rows or columns exchange tokens during
  randomization. Supported for all sequential algorithms (curvecat, swapcat,
  tswapcat) and flows through to quantize(), batch functions, trace
  diagnostics, and vegan integration.

# nullcat 0.1.0

* Initial CRAN submission.
* Categorical generalizations of binary null models (curvecat, swapcat, tswapcat, r0cat, c0cat)
* Stratified quantitative null models via quantize()
* MCMC diagnostics (trace_cat, suggest_n_iter)
* Integration with vegan package

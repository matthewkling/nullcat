# nullcat (development version)

* New `row_groups` and `col_groups` arguments randomize within groups of rows and/or columns. The matrix is split into blocks (every combination of a row group and a column group) and each block is randomized independently, so values never leave their block. For fixed-margin algorithms this gives a constrained null model, e.g. one in which each region keeps its own species pool. Supported by every algorithm, `nullcat()`, `quantize()`, `curvecat()`, `swapcat()`, `tswapcat()`, the batch functions, `trace_cat()`, `suggest_n_iter()`, and the `vegan` commsim functions. In `quantize()`, strata are defined once for the whole matrix and quantitative values are reassigned within each block. When groups are used, `n_iter` is applied to every block, and `trace_cat()` / `suggest_n_iter()` trace the largest block.

* `nullcat_batch()` and `quantize_batch()` also accept a matrix of groups with one column per replicate, using column `r` for replicate `r`.

* New function `tile_groups()` generates random hexagonal tilings of sites, one per column, for spatially constrained null models: passed to the batch functions, each replicate is randomized within the tiles of a different, randomly placed tiling, so sites are pooled with their neighbours. With longitude/latitude coordinates, tiles are laid out on the sphere (a randomly rotated Fibonacci lattice) with `spacing` in metres; planar coordinates use a randomly rotated and offset hexagonal grid. Two sites share a tile in about half of the tilings at a distance of `0.42 * spacing`.

* New function `group_sizes()` returns the average size of each row's group across a grouping vector or matrix, to check pool sizes before randomizing (e.g. while choosing a tile `spacing`).

* `wt_row` and `wt_col` are deprecated. Weighting which pairs of rows or columns exchange values does not change the null distribution that `curvecat` and `swapcat` converge to, only how quickly they mix: with weights that connect all rows (e.g. distance-decay kernels), results converge to the unconstrained null as `n_iter` grows. The only lasting effect of a weight matrix is to separate rows that are not connected by any chain of nonzero weights, so weight matrices are now converted to `row_groups` / `col_groups` (the connected components of their nonzero entries) with a warning; if the weights connect everything they have no effect, and if they have no off-diagonal entries an error is raised. This also fixes a bug: under `quantize(fixed = "cell")`, row weights were ignored on half the iterations, so values could move between supposedly separated groups.

* `quantize_prep()` (and so `quantize()` and its relatives) now warns when `fixed = "cell"` and some rows or columns have all their values in a single stratum. Such rows and columns can never change under cell mode, because cells only move by trading places with cells of a different stratum. Using `zero_stratum = TRUE` or a `transform` usually separates low values into their own strata.

* Bug fix: with a single row or column, or `n_iter = 0`, the sequential algorithms (`curvecat`, `swapcat`, `tswapcat`) returned the input categories even when `output = "index"` was requested, so `quantize(fixed = "cell")` scrambled values for such matrices. They now return the identity index.

* Documentation corrections: `curvecat` and `swapcat` converge to the same (uniform) null distribution, but `tswapcat` does not; its redraw-until-valid design biases sampling toward matrices with many available swaps (Miklós & Podani 2004).

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

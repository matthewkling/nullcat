# nullcat

`nullcat` provides null model algorithms for categorical and
quantitative community ecology data. It extends classic binary null
models (e.g., curveball, swap) to work with categorical data, and
introduces a stratified randomization framework for continuous data that
addresses limitations in existing methods for randomizing quantitative
and count data.

## Installation

``` r
# Install stable version from CRAN:
install.packages("nullcat")

# Or dev version from GitHub:
# install.packages("remotes")
remotes::install_github("matthewkling/nullcat")
```

## Quick start

#### Categorical null models

Generalize binary null models to matrices where cells contain integers
representing discrete categories:

``` r
library(nullcat)

# Create a categorical matrix
set.seed(123)
cat_matrix <- matrix(sample(1:4, 20*10, replace = TRUE), nrow = 20)

# Randomize using curvecat (preserves row & column category multisets)
randomized <- curvecat(cat_matrix, n_iter = 1000)

# Verify margins are preserved
all.equal(sort(cat_matrix[1,]), sort(randomized[1,]))
#> [1] TRUE
all.equal(sort(cat_matrix[,1]), sort(randomized[,1]))
#> [1] TRUE
```

Available algorithms:
[`curvecat()`](https://matthewkling.github.io/nullcat/reference/curvecat.md),
[`swapcat()`](https://matthewkling.github.io/nullcat/reference/swapcat.md),
[`tswapcat()`](https://matthewkling.github.io/nullcat/reference/tswapcat.md),
[`r0cat()`](https://matthewkling.github.io/nullcat/reference/r0cat.md),
[`c0cat()`](https://matthewkling.github.io/nullcat/reference/c0cat.md)

#### Quantitative null models

Apply stratified randomization to continuous community data:

``` r
# Create a quantitative community matrix
set.seed(456)
comm <- matrix(rexp(50 * 30, rate = 0.5), nrow = 50)

# Stratified randomization with 5 strata
rand1 <- quantize(comm, n_strata = 5, n_iter = 2000)

# Preserve row value multisets (row sums maintained)
rand2 <- quantize(comm, n_strata = 5, fixed = "row", n_iter = 2000)
all.equal(rowSums(comm), rowSums(rand2))
#> [1] TRUE
```

## Learn more

See
[`vignette("nullcat")`](https://matthewkling.github.io/nullcat/articles/nullcat.md)
for comprehensive documentation including:

- Categorical null model algorithms and their constraints
- Quantitative null model workflow and stratification options
- Convergence diagnostics and burn-in estimation
- Efficient batch generation of null distributions
- Integration with the vegan package

**Package website:** <https://matthewkling.github.io/nullcat/>

**Report issues:** <https://github.com/matthewkling/nullcat/issues>

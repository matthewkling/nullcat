## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(nullcat)

## ----categorical-basic--------------------------------------------------------
# Create a categorical matrix
set.seed(123)
x <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)

# Randomize using curvecat (preserves row & column category multisets)
x_rand <- curvecat(x, n_iter = 1000)

# Verify margins are preserved
all.equal(sort(x[1, ]), sort(x_rand[1, ]))  # Row preserved
all.equal(sort(x[, 1]), sort(x_rand[, 1]))  # Column preserved

## ----categorical-margins------------------------------------------------------
set.seed(456)
x <- matrix(sample(1:3, 60, replace = TRUE), nrow = 10)

# Fixed-fixed: both row and column margins preserved
x_fixed <- curvecat(x, n_iter = 1000)

# Row-constrained: only row margins preserved
x_row <- r0cat(x)

# Column-constrained: only column margins preserved
x_col <- c0cat(x)

# Check row sums by category
table(x[1, ])           # Original row 1
table(x_fixed[1, ])     # Preserved in curvecat
table(x_row[1, ])       # Preserved in r0cat
table(x_col[1, ])       # Not preserved in c0cat

## ----quantitative-basic-------------------------------------------------------
# Create a quantitative community matrix
set.seed(789)
comm <- matrix(runif(200), nrow = 20)

# Default: curvecat-backed stratified randomization with 5 strata
rand1 <- quantize(comm, n_strata = 5, n_iter = 2000, fixed = "row")

# Values are similar but rearranged
cor(as.vector(comm), as.vector(rand1))
plot(rowSums(comm), rowSums(rand1))
plot(colSums(comm), colSums(rand1))



## ----quantitative-strata------------------------------------------------------
set.seed(200)
x <- rexp(100, .1)

# More strata = less mixing, higher fidelity to original distribution
s3 <- stratify(x, n_strata = 3)
s10 <- stratify(x, n_strata = 10)

table(s3)   # Coarser bins
table(s10)  # Finer bins

# Transform before stratifying (e.g., log-transform for skewed data)
s_log <- stratify(x, n_strata = 5, transform = log1p)
table(s_log)

# Rank transform creates equal-occupancy strata
s_rank <- stratify(x, n_strata = 5, transform = rank)
table(s_rank)  # Nearly equal counts per stratum

# Separate zeros into their own stratum
x_with_zeros <- c(0, 0, 0, x)
s_zero <- stratify(x_with_zeros, n_strata = 4, zero_stratum = TRUE)
table(s_zero)

## ----quantitative-fixed-------------------------------------------------------
set.seed(100)
comm <- matrix(rexp(100), nrow = 10)

# Preserve row value multisets (quantitative row sums maintained)
rand_row <- quantize(comm, n_strata = 5, fixed = "row", n_iter = 2000)
all.equal(rowSums(comm), rowSums(rand_row))

# Preserve column value multisets (quantitative column sums maintained)
rand_col <- quantize(comm, n_strata = 5, fixed = "col", n_iter = 2000)
all.equal(colSums(comm), colSums(rand_col))

# Cell-level preservation: each value moves with its original cell location
# The categorical randomization determines WHERE cells go, but each cell
# carries its original value with it
rand_cell <- quantize(comm, n_strata = 5, fixed = "cell", n_iter = 2000)

# Values shuffled globally within strata, holding none of the above fixed
rand <- quantize(comm, n_strata = 5, fixed = "stratum", n_iter = 2000)

# For non-fixed-fixed methods like r0cat or c0cat, only some fixed options make sense:
# r0cat with fixed="col" or c0cat with fixed="row" would be incompatible

## ----efficient-batch----------------------------------------------------------
set.seed(400)
comm <- matrix(rexp(200), nrow = 20)

# Prepare once
prep <- quantize_prep(comm, n_strata = 5, fixed = "row", n_iter = 2000)

# Generate many randomizations efficiently
rand1 <- quantize(prep = prep)
rand2 <- quantize(prep = prep)
rand3 <- quantize(prep = prep)

# Or use batch functions
nulls <- quantize_batch(comm, n_reps = 99, n_strata = 5, 
                       fixed = "row", n_iter = 2000)
dim(nulls)  # 20 rows × 10 cols × 99 replicates

## ----trace-diagnostics--------------------------------------------------------
set.seed(300)
x <- matrix(sample(1:5, 400, replace = TRUE), nrow = 20)

# Generate trace showing mixing over iterations
trace <- trace_cat(x, fun = "nullcat", method = "curvecat",
                   n_iter = 1000, n_chains = 3, thin = 10)

# Visual inspection
plot(trace)

# Automatic burn-in suggestion
suggested <- suggest_n_iter(trace, tail_frac = 0.3)
print(suggested)

## ----vegan-integration, eval = FALSE------------------------------------------
# library(vegan)
# 
# # Categorical data
# x_cat <- matrix(sample(1:4, 100, replace = TRUE), nrow = 10)
# cs <- nullcat_commsim_seq(method = "curvecat")
# nm <- nullmodel(x_cat, cs)
# sims <- simulate(nm, nsim = 99, burnin = 1000, thin = 100)
# 
# # Quantitative data
# x_quant <- matrix(rexp(100), nrow = 10)
# cs_quant <- quantize_commsim(n_strata = 5, method = "curvecat", n_iter = 2000)
# nm_quant <- nullmodel(x_quant, cs_quant)
# sims_quant <- simulate(nm_quant, nsim = 99)


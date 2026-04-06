pkgname <- "phylospatial"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('phylospatial')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("clade_dist")
### * clade_dist

flush(stderr()); flush(stdout())

### Name: clade_dist
### Title: Pairwise distances among clades or nodes
### Aliases: clade_dist

### ** Examples

if(requireNamespace("phytools", quietly = TRUE)){
  clade_dist(ape::rtree(10))
}



cleanEx()
nameEx("moss")
### * moss

flush(stderr()); flush(stdout())

### Name: moss
### Title: Load California moss spatial phylogenetic data
### Aliases: moss

### ** Examples





cleanEx()
nameEx("phylospatial")
### * phylospatial

flush(stderr()); flush(stdout())

### Name: phylospatial
### Title: Create a spatial phylogenetic object
### Aliases: phylospatial

### ** Examples




cleanEx()
nameEx("plot.phylospatial")
### * plot.phylospatial

flush(stderr()); flush(stdout())

### Name: plot.phylospatial
### Title: Plot a 'phylospatial' object
### Aliases: plot.phylospatial

### ** Examples

ps <- ps_simulate(20, 20, 20)
plot(ps, "tree")
plot(ps, "comm")



cleanEx()
nameEx("plot_lambda")
### * plot_lambda

flush(stderr()); flush(stdout())

### Name: plot_lambda
### Title: Plot alternative lambda values
### Aliases: plot_lambda

### ** Examples

plot_lambda()
plot_lambda(seq(0, 3, .1))




cleanEx()
nameEx("ps_add_dissim")
### * ps_add_dissim

flush(stderr()); flush(stdout())

### Name: ps_add_dissim
### Title: Add community dissimilarity data to a 'phylospatial' object
### Aliases: ps_add_dissim

### ** Examples

ps <- ps_simulate(data_type = "prob")
ps_add_dissim(ps)
ps_add_dissim(ps, fun = "vegdist", method = "jaccard", endemism = TRUE)




cleanEx()
nameEx("ps_canape")
### * ps_canape

flush(stderr()); flush(stdout())

### Name: ps_canape
### Title: Categorical Analysis of Neo- and Paleo-Endemism (CANAPE)
### Aliases: ps_canape

### ** Examples





cleanEx()
nameEx("ps_canaper")
### * ps_canaper

flush(stderr()); flush(stdout())

### Name: ps_canaper
### Title: Binary randomization tests including CANAPE
### Aliases: ps_canaper

### ** Examples




cleanEx()
nameEx("ps_dissim")
### * ps_dissim

flush(stderr()); flush(stdout())

### Name: ps_dissim
### Title: Quantitative phylogenetic dissimilarity
### Aliases: ps_dissim

### ** Examples

# example data set:
ps <- ps_simulate(n_tips = 50)

# The default arguments give Sorensen's quantitative dissimilarity index
# (a.k.a. Bray-Curtis distance):
d <- ps_dissim(ps)

# Specifying a custom formula explicitly via `designdist`;
# (this is the Bray-Curtis formula, so it's equivalent to the prior example)
d <- ps_dissim(ps, method = "(b+c)/(2*a+b+c)",
      fun = "designdist", terms = "minimum", abcd = TRUE)

# Alternative arguments can specify a wide range of dissimilarity measures;
# here's endemism-weighted Jaccard's dissimilarity:
d <- ps_dissim(ps, method = "jaccard", endemism = TRUE)




cleanEx()
nameEx("ps_diversity")
### * ps_diversity

flush(stderr()); flush(stdout())

### Name: ps_diversity
### Title: Calculate spatial phylogenetic diversity and endemism metrics
### Aliases: ps_diversity

### ** Examples

ps <- ps_simulate()
div <- ps_diversity(ps)
terra::plot(div)




cleanEx()
nameEx("ps_expand")
### * ps_expand

flush(stderr()); flush(stdout())

### Name: ps_expand
### Title: Expand occupied-only results to full spatial extent
### Aliases: ps_expand

### ** Examples

ps <- ps_simulate()

# custom analysis on the occupied-only community matrix
site_totals <- matrix(rowSums(ps$comm), ncol = 1)
colnames(site_totals) <- "total"

# expand to full extent as a matrix
ps_expand(ps, site_totals)

# expand and convert to spatial
ps_expand(ps, site_totals, spatial = TRUE)



cleanEx()
nameEx("ps_geodist")
### * ps_geodist

flush(stderr()); flush(stdout())

### Name: ps_geodist
### Title: Geographic distance between sites
### Aliases: ps_geodist

### ** Examples




cleanEx()
nameEx("ps_get_comm")
### * ps_get_comm

flush(stderr()); flush(stdout())

### Name: ps_get_comm
### Title: Get 'phylospatial' community data
### Aliases: ps_get_comm

### ** Examples

ps <- ps_simulate()

# the defaults return a spatial object of terminal taxa distributions:
ps_get_comm(ps)

# get distributions for all taxa, as a matrix
pcomm <- ps_get_comm(ps, tips_only = FALSE, spatial = FALSE)




cleanEx()
nameEx("ps_ordinate")
### * ps_ordinate

flush(stderr()); flush(stdout())

### Name: ps_ordinate
### Title: Community phylogenetic ordination
### Aliases: ps_ordinate

### ** Examples

ps <- ps_add_dissim(ps_simulate(50, 5, 5))
ord <- ps_ordinate(ps, method = "cmds", k = 4)
terra::plot(ord)




cleanEx()
nameEx("ps_prioritize")
### * ps_prioritize

flush(stderr()); flush(stdout())

### Name: ps_prioritize
### Title: Phylogenetic conservation prioritization
### Aliases: ps_prioritize

### ** Examples




cleanEx()
nameEx("ps_quantize")
### * ps_quantize

flush(stderr()); flush(stdout())

### Name: ps_quantize
### Title: Stratified randomization of a phylospatial object
### Aliases: ps_quantize

### ** Examples




cleanEx()
nameEx("ps_rand")
### * ps_rand

flush(stderr()); flush(stdout())

### Name: ps_rand
### Title: Null model randomization analysis of alpha diversity metrics
### Aliases: ps_rand

### ** Examples




cleanEx()
nameEx("ps_regions")
### * ps_regions

flush(stderr()); flush(stdout())

### Name: ps_regions
### Title: Cluster analysis to identify phylogenetic regions
### Aliases: ps_regions

### ** Examples

ps <- ps_simulate()

# using kmeans clustering algorithm
terra::plot(ps_regions(ps, method = "kmeans"))

# to use a hierarchical clustering method, first we have to `ps_add_dissim()`
terra::plot(ps_regions(ps_add_dissim(ps), k = 7, method = "average"))




cleanEx()
nameEx("ps_regions_eval")
### * ps_regions_eval

flush(stderr()); flush(stdout())

### Name: ps_regions_eval
### Title: Evaluate region numbers
### Aliases: ps_regions_eval

### ** Examples

ps <- ps_add_dissim(ps_simulate())
ps_regions_eval(ps, k = 1:15, plot = TRUE)




cleanEx()
nameEx("ps_rgb")
### * ps_rgb

flush(stderr()); flush(stdout())

### Name: ps_rgb
### Title: Map phylospatial data onto RGB color bands
### Aliases: ps_rgb

### ** Examples

ps <- ps_add_dissim(ps_simulate(50, 20, 20))
RGB <- ps_rgb(ps, method = "cmds")
terra::plotRGB(RGB * 255, smooth = FALSE)




cleanEx()
nameEx("ps_simulate")
### * ps_simulate

flush(stderr()); flush(stdout())

### Name: ps_simulate
### Title: Simulate a toy spatial phylogenetic data set
### Aliases: ps_simulate

### ** Examples

# using all the defaults
ps_simulate()

# specifying some arguments
plot(ps_simulate(n_tips = 50, n_x = 30, n_y = 40, data_type = "abundance"), "comm")



cleanEx()
nameEx("quantize")
### * quantize

flush(stderr()); flush(stdout())

### Name: quantize
### Title: Stratified randomization of community matrix
### Aliases: quantize

### ** Examples




cleanEx()
nameEx("to_spatial")
### * to_spatial

flush(stderr()); flush(stdout())

### Name: to_spatial
### Title: Convert a site-by-variable matrix into a SpatRaster or sf object
### Aliases: to_spatial

### ** Examples

ps <- moss()
# ps$comm contains only occupied sites, so expand before converting:
to_spatial(ps_expand(ps, ps$comm[, 1:5]), ps$spatial)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

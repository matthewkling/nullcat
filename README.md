
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nullcat

Null models for categorical and quantized continuous data in community
ecology.

This package extends classic **binary** null models to work with
**categorical** data. Categorical generalizations currently include
`curveball -> curvecat` and `swap -> swapcat`.

It also provides a routine for using these algorithms with
**continuous** data. This `quantize()` routine works by converting
continuous values to discrete strata, randomizing them using categorical
algorithms, and mapping the randomized categories back to the original
quantities.

# Steps to reproduce simulations and networks "Graphical Models for Zero-Inflated Single Cell Gene Expression"

Install `HurdleNormal` (eg with `devtools::install_github`) with
`dependencies = TRUE` and install `netbenchmark`, `Mus.musculus` from bioconductor.

For simulations:
1. source simulations.R
2. Run processSimulations.Rmd

For networks generated from Shalek et al data:

1. Run fitAlex.R
2. Run processShalekMAITNetworks.Rmd

For networks generated from TFh data, the source data is as-yet unpublished a primary journal so cannot yet be made available publically.

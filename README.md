# Digital sampling and distance data

This repository contains an R simulation study comparing population-detection estimates for digital surveys and traditional distance sampling. The article will describe the framework: detection can vary with distance, covariates, and their interaction, and assuming perfect detection can bias estimates.

## Workflow

1. Run a `Dataset_creation_regime*` script to simulate survey data under one of three distance-effect regimes: half-normal, beta-shaped, or uniform. The scripts save detected-animal data and simulation summaries.
2. Run `fit_all_approaches_optimizing_C.R` to fit the estimation approaches to those data and write results to CSV. It compares distance-only, covariate, stratified, beta/normal, uniform-detection, and true-shape estimates.

More detailed descriptions to come. 



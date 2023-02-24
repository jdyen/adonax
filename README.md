## Analysis of the effects of *Arundo donax* on fish assemblages in north-eastern Spain

This repository contains code and data to support the following manuscript:
Maceda-Viega A, Mac Nally R, Cano-Rocabayera O, de Sostoa A, and Yen JDL. Giant reed (*Arundo donax*) invasion may have little effects on fish assemblages in multi-stressed streams. In preparation.


Last updated: 24 February 2023 


### Overview

This study used an extensive data set spanning north-eastern Spain to address two aims: (1) identify the environmental characteristics of the stream reaches invaded by *A. donax* in north-eastern Spain to contextualize potential impacts on fishes; and (2) to establish associations between *A. donax* and the relative densities, taxonomic richness, and individual body condition of fishes. We used hierarchical Bayesian models to estimate associations with environmental predictor variables and *A. donax*. This repository contains all code to reproduce these analyses and associated outputs included in the manuscript.
 

### Usage

The `main.R` script runs the primary analysis and generates associated outputs. Several helper functions are provided in the `utils.R` script in the `R` directory; these are sourced directly within `main.R`. All models are fitted with the Stan software called via the `rstan` R package. Stan models are in the `src` directory. Several additional directories will need to be created to store model outputs: `outputs/figures/`, `outputs/fitted/`, and `outputs/tables/`.

Data to reproduce analyses are not currently included in this repository but will be uploaded following manuscript submission.


### Contact

For additional information about the analyses in this repository, please contact Jian Yen (jian.yen [at] delwp.vic.gov.au).


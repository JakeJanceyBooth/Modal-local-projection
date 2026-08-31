# Modal Local Projections Monte Carlo Code

## Overview

Linear modal local projection can be viewed as a scalar, single-equation specialization of the 
dynamic vector mode regression of Kemp et al. (2020). We use simulations to examine its finite-sample 
performance in forecasting and dynamic response analysis, along with blocks-of-blocks bootstrap 
inference, and compare the results to those from mean- and median-based methods.

## Repository structure

- `R/` contains reusable estimators and DGPs
- `scripts/` contains the Monte Carlo experiments and figure creation
- `results/` contains selected simulation output/summaries and final figures
- `demo.R` is a quick example of the estimators and DGPs

## Requirements

The analysis uses R and the following R packages:

- `data.table`
- `dplyr`
- `ggplot2`
- `quantreg`
- `tidyr`
- `vars`

## Running the code

All scripts assume that the working directory is the repository root. Paths are
constructed relative to that directory using `file.path()`.

Final run-settings objects record the R and package environment where
appropriate using `sessionInfo()`.

## Experiments

The scripts are by default set to computationally-inexpensive
settings, in particular only R = 10 Monte Carlo replications. At the number of
replications used in the paper (e.g. R = 1000, 5000), the scripts can take several 
hours to complete a full run.

Monte Carlo summaries are calculated over successful model fits or forecasts.
Success counts and rates are reported separately. Failed numerical fits are 
treated as unsuccessful observations and excluded .

## Results and figures

The `results/` folder contains selected outputs from the Monte Carlo
simulations.

- Final summary objects used for analysis and figures are saved.
- Certain objects are saved to allow later summaries/statistics to be
  produced without having to rerun expensive simulations.
- Large raw forecast tables and exploratory outputs are not saved.
- The final figure PDFs are produced in 05_make_figures.R and saved 
  in `results/figures/`.














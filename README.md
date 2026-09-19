# Modal Local Projections Monte Carlo Code

## Overview

Linear modal local projection can be viewed as a scalar, single-equation specialization of the 
dynamic vector mode regression of Kemp et al. (2020). We use Monte Carlo simulations to examine its finite-sample 
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

## Experiments

Scripts 01–04 run the Monte Carlo experiments with 10 replications by default. 
The saved results used in the paper correspond to 100, 1,000, 5,000, and 100 
replications, respectively. These settings can be restored by changing n_replications 
in scripts 01–03 and n_outer_replications in script 04. Full runs could take several hours.

Monte Carlo summaries are calculated over successful model fits or forecasts.
Success counts and rates are reported separately. Failed numerical fits are 
treated as unsuccessful observations and excluded.

## Results and figures

The `results/` folder contains selected outputs from the Monte Carlo
simulations.

- Final summary objects used for analysis and figures are saved.
- Certain objects are saved to allow later summaries/statistics to be
  produced without having to rerun expensive simulations.
- Large raw forecast tables and exploratory outputs are not saved.
- Running scripts/05_make_figures.R recreates the figure PDFs in 
  results/figures/ using the saved simulation outputs already included 
  in the repository, without rerunning any simulations.














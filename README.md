# Enhancing Regional Sea Level Predictions: A Data-Driven Dual-Module Approach

This repository contains the code and data processing scripts for the research project.

## Environment Setup

This project requires both Python and R environments. Follow the instructions below to set up the environment.

### Python Setup

It is recommended to use `conda` to manage the Python environment.

1. Create a new conda environment:

   ```bash
   conda create -n sea_level_env python=3.11
   conda activate sea_level_env
   ```

2. Install required packages using the provided setup script:

   ```bash
   bash setup.sh
   ```

### R Setup

1. Install R (version 4.5 or higher).
2. Install required R packages using the provided requirements script:

   ```R
   source("requirements.R")
   ```

## Data Preparation

1. Download ERA5 monthly data:

   ```bash
   python utils/download_era5_monthly.py
   ```

2. Download IPCC AR6 sea level and VLM projections:

   ```bash
   python utils/ar6.py
   ```

Additional data files will be automatically downloaded when required by the scripts.

## Results Reproduction

1. To reproduce the main results from the paper:

   ```R
   source("update_results.R")
   ```

   ```bash
   bash generate/generate_cluster_prediction_diff.sh
   ```

2. To generate the complete backtesting and prediction results for supplementary materials:

   ```R
   source("generate_supplementary.R")
   ```

All results will be saved in the `results/` directory.
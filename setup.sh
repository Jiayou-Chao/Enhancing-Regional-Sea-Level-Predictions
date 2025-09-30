#!/bin/bash

conda install -c conda-forge xarray dask netCDF4 bottleneck requests -y
conda install -c conda-forge eccodes -y
conda install -c conda-forge tqdm -y
conda install scipy -y
pip install cfgrib
conda install hdbscan plotly -y
pip install kaleido
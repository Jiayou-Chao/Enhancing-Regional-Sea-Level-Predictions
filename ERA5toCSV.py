import argparse
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


def collect_annual_averages(station: str) -> None:
    station_dir = Path('C:/Users/chaoj/OneDrive/R/RSL_prediction/dataset/era5') / station
    if not station_dir.exists():
        raise FileNotFoundError(f"ERA5 directory not found for station '{station}'")

    for subdir in sorted(d for d in station_dir.iterdir() if d.is_dir()):
        files = sorted(f for f in subdir.iterdir() if f.is_file())

        annual_averages = {'year': [], subdir.name: []}

        for file in files:
            year = file.name[:4]
            if year == '2024':
                continue

            with nc.Dataset(file) as dataset:
                variable = list(dataset.variables.keys())[-1]
                annual_averages['year'].append(year)
                annual_averages[subdir.name].append(np.mean(dataset.variables[variable][:]))

        output_dir = Path('data/plots') / station / 'era5data'
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f'{subdir.name}_average.csv'
        pd.DataFrame(annual_averages).to_csv(output_file, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description='Convert ERA5 station data to annual averages.')
    parser.add_argument('station', type=str, help='Tide gauge station name.')
    args = parser.parse_args()

    collect_annual_averages(args.station)


if __name__ == '__main__':
    main()

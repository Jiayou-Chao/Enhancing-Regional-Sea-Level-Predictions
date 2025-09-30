import xarray as xr
import pandas as pd
import os
from pathlib import Path
from tqdm import tqdm
import sys
import pathlib
import argparse
import numpy as np
import multiprocessing

sys.path.append(os.path.join(pathlib.Path(__file__).parent.resolve()))
import ar6

era_dir = "data/era5data"

def grib_to_df(era_dir: str, var_name: str, longitude_of_interest: float, latitude_of_interest: float, verbose: bool = False, save: bool = False, save_dir: str = '') -> pd.DataFrame:
    data_path = os.path.join(era_dir, var_name)

    all_files = os.listdir(data_path)
    grib_files = [file for file in all_files if file.endswith(".grib")]
    avg_dict = {"Year": [], var_name: []}
    for i in tqdm(range(len(grib_files))):
        file_path = os.path.join(data_path, grib_files[i])
        ds = xr.open_dataset(file_path, engine='cfgrib')
        var = list(ds.data_vars)[0]
        year = int(ds.time.dt.year.values[0])
        ds = ds.mean(dim = "time")
        avg = float(ds[var].sel(longitude = longitude_of_interest,latitude=latitude_of_interest,method = "bfill").values)
        avg_dict["Year"].append(year)

        if np.isnan(avg):
            non_nan_mask = ~np.isnan(ds[var])

            lats = ds.latitude.values
            lons = ds.longitude.values

            lons, lats = np.meshgrid(lons, lats)

            distances = np.sqrt((lats - latitude_of_interest)**2 + (lons - longitude_of_interest)**2)
            distances = np.where(non_nan_mask, distances, np.inf)

            min_dist_index = np.unravel_index(np.argmin(distances), distances.shape)


            avg = float(ds[var].values[min_dist_index])



        avg_dict[var_name].append(avg)
        if verbose:
            print(f"Year {year} has an average of {avg}")

    avg_df = pd.DataFrame(avg_dict)

    if save:
        avg_df.to_csv(os.path.join(save_dir, f"{var_name}_avg.csv"), header=True)

    for idx_file in Path(data_path).glob('*.idx'):
        idx_file.unlink()

    return avg_df

def process_var(args):
    var, era_dir, longitude_of_interest, latitude_of_interest = args
    print("Starting " + var)
    return grib_to_df(era_dir, var, longitude_of_interest, latitude_of_interest)

def get_era5_from_coords(longitude_of_interest: float, latitude_of_interest: float, era_dir: str = "data/era5data", save: bool = True, save_path: str = None) -> pd.DataFrame:
    variable_names = os.listdir(era_dir)
    out = pd.DataFrame({"Year": list(range(1950, 2024))})

    with multiprocessing.Pool() as pool:
        results = pool.map(process_var, [(var, era_dir, longitude_of_interest, latitude_of_interest) for var in variable_names])

    for var_df in results:
        out = pd.merge(out, var_df, on="Year")

    if save and save_path is not None:
        save_dir = os.path.dirname(save_path)
        os.makedirs(save_dir, exist_ok=True)
        out.to_csv(save_path, header=True, index=False)

    return out

def get_era5_from_station(station: int | str, save_dir: str = 'data/era5data_processed', from_file: bool = True) -> pd.DataFrame:
    station_id = ar6.get_station_id(station)
    latlon = ar6.station2lonlat(station_id)
    lat = latlon['lat']
    lon = latlon['lon']
    if save_dir is not None:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        save = True
        save_path = os.path.join(save_dir, f"{str(station_id)}_era5_data.csv")
    if from_file and os.path.exists(save_path):
        return pd.read_csv(save_path)
    return get_era5_from_coords(longitude_of_interest=lon, latitude_of_interest=lat, save=save, save_path=save_path)

def save_era5_stations(stations: list[int], **kwargs) -> None:
    for station in stations:
        get_era5_from_station(station, from_file=False, **kwargs)

def delete_idx_files(directory: str = "data/era5data") -> None:
    """Delete all .idx files in the specified directory."""
    for idx_file in Path(directory).rglob('*.idx'):
        idx_file.unlink()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process ERA5 data for given station IDs.')
    parser.add_argument('stations', nargs='+', type=int, help='Station IDs to process')
    parser.add_argument('-u', '--update', action='store_true', help='Force update the data file')
    args = parser.parse_args()

    for station in args.stations:
        get_era5_from_station(station, from_file=not args.update)

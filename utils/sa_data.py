import xarray as xr
import pandas as pd
import os
import glob
from tqdm import tqdm
from datetime import datetime
import numpy as np

import sys
import pathlib
sys.path.append(os.path.join(pathlib.Path(__file__).parent.resolve()))

import ar6

def list_files(folder_path, extension='*', recursive=False, full_path=True, **kwargs):
    if extension == '*':
        pattern = '**/*' if recursive else '*'
    else:
        pattern = f'**/*.{extension}' if recursive else f'*.{extension}'

    files = glob.glob(os.path.join(folder_path, pattern), recursive=recursive, **kwargs)
    files.sort()

    if full_path:
        return [os.path.abspath(file) for file in files]
    return files



def list_satellite_files(**kwargs):
    return list_files('data/satellite/monthly_raw', 'nc', recursive=False, **kwargs)

def read_one_satellite_data(file_path, lat, lon):
    """
    It will search for the lon_bnds and lat_bnds in the file to locate the data point.
    The data files are from https://cds.climate.copernicus.eu/cdsapp
    data['lon_bnds'].shape (1440, 2)
    data['lat_bnds'].shape (720, 2)
    After locating the lon and lat, it will return the sla variable.
    data['sla'].shape data['sla'].shape
    """
    data = xr.open_dataset(file_path)
    lon_bnds = data['lon_bnds']
    lat_bnds = data['lat_bnds']
    if lon < 0:
        lon += 360
    lon_idx = np.where((lon >= lon_bnds[:, 0]) & (lon <= lon_bnds[:, 1]))[0][0]
    lat_idx = np.where((lat <= lat_bnds[:, 0]) & (lat >= lat_bnds[:, 1]))[0][0]
    sla = data['sla'][0, lat_idx, lon_idx]
    return sla

def read_satellite_data(lat, lon, return_dataframe=True,
                        save_path=None):
    files = list_satellite_files()
    dates = []
    data = []
    for file in tqdm(files):
        sla = read_one_satellite_data(file, lat=lat, lon=lon)
        data.append(sla.item())
        dates.append(pd.to_datetime(sla.time.values).strftime('%Y-%m-%d'))
    result = {'date': dates, 'sla': data}
    if return_dataframe:
        result = pd.DataFrame(result)
        if 'sla' not in result.columns:
            result['sla'] = np.nan
        if save_path:
            result.to_csv(save_path, index=False)
    return result

def read_satellite_data_station(station_name: int|str, save_dir="data/satellite/processed", from_file: bool = True, **kwargs):
    try:
        station_id = int(station_name)
    except ValueError:
        if station_name in ar6.station_loc_id_map:
            station_id = ar6.station_loc_id_map[station_name]
        else:
            station_id = ar6.get_ar6_station_id(station_name)

    save_path = os.path.join(save_dir, f"{station_id}.csv")
    os.makedirs(save_dir, exist_ok=True)

    if from_file and os.path.exists(save_path):
        return pd.read_csv(save_path)
    latlon = ar6.station2lonlat(station_id)
    lat, lon = latlon['lat'], latlon['lon']
    return read_satellite_data(lat=lat, lon=lon, save_path=save_path, **kwargs)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.description = "This script reads and processes satellite data and save them to 'data/satellite/processed'."
    parser.add_argument('station_ids', type=str, nargs='+', help='One or more station IDs for sea level analysis')
    parser.add_argument('-u', '--update', action='store_true', help='Update the files.')
    args = parser.parse_args()
    from_file = not args.update
    for station_id in args.station_ids:
        read_satellite_data_station(station_id, from_file=from_file)

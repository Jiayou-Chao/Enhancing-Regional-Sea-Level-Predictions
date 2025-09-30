import os
import sys

import cdsapi
from tqdm import tqdm

from data import station2lonlat

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dotenv import load_dotenv

import utils

load_dotenv()

def download_era5(lat, lon, year, var: str|list, save_path, client=cdsapi.Client()):

    client.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': var,
            'year': year,
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': [
                lat + 0.001, lon - 0.001, lat - 0.001, lon + 0.001,
            ],
        },
        save_path)

def download_era5_station(station_name, year, var, save_path, **kwargs):
    if station_name not in utils.possible_station_names:
        print("Possible station names:")
        print(utils.possible_station_names)
        station_name = input("Enter a valid station name: ")
    latlon = station2lonlat(station_name)
    lat, lon = latlon['lat'], latlon['lon']
    download_era5(lat, lon, year, var, save_path)

def main():
    dataset_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'dataset/era5'))
    remote_folder = "gds:RSL_prediction/dataset/era5"
    years = range(1940, 2025)
    variables = ['2m_temperature', 'mean_total_precipitation_rate', '2m_dewpoint_temperature', 'mean_sea_level_pressure']
    station = 'monterey'
    for year in tqdm(years):
        for var in variables:
            save_dir = os.path.join(dataset_folder, station, var)
            save_path = os.path.join(save_dir, f'{year}.nc')
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            if not os.path.exists(save_path):
                print(f"Downloading {var} for {station} in {year}...")
                download_era5_station(station, year, var, save_path)
    os.system(f"rclone copy {dataset_folder} {remote_folder} --size-only -u --checkers 16")


if __name__ == '__main__':
    main()

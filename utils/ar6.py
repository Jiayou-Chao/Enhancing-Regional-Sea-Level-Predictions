
import argparse
import xarray as xr
import os
import pandas as pd
import requests
import numpy as np
from scipy.interpolate import interp1d

station_loc_id_map = {
    'astoria': 256,
    'battery_park': 12,
    'crescent_city': 378,
    'south_beach': 1196,
    'monterey': 1352,
    'newport': 351,
    'tofino': 165,
    'cape_charles': 636,
    'charleston': 234,
    'elly_oil_platform': 245
}

global_station_latlon = {
    "london": (51.5074, -0.1278),
    "tokyo": (35.6895, 139.6917),
    "beijing": (39.9042, 116.4074),
    "new_york": (40.7128, -74.0060),
    "sydney": (-33.8688, 151.2153)
}


def load_ar6_data(scenario, workflow_ids=['wf_1e', 'wf_1f', 'wf_2e', 'wf_2f', 'wf_3e', 'wf_3f', 'wf_4']):
    """
    Load AR6 sea level projection data for given scenarios and workflow IDs.

    Args:
    scenario (str): The SSP scenario to load.
    workflow_ids (list): The list of workflow IDs to load. Read more about it at https://github.com/Rutgers-ESSP/IPCC-AR6-Sea-Level-Projections/blob/main/FAQ.md

    Returns:
    dict: A dictionary of xarray.Datasets keyed by workflow ID.
    """
    datasets = {}
    for workflow_id in workflow_ids:
        base_url = f'https://storage.googleapis.com/ar6-lsl-simulations-public-standard/tide-gauges/full_sample_workflows/{workflow_id}'
        data_url = f'{base_url}/{scenario}/total-workflow.zarr'
        datasets[workflow_id] = xr.open_dataset(data_url, engine='zarr', chunks=None)
    return datasets

def get_median_sea_level_change(ds, location_id, use_2020_as_reference=True):
    """
    Calculate median sea level change for a specific location relative to 2020.

    Args:
    ds (xarray.Dataset): The dataset containing sea level projections.
    location_id (int): The ID of the location to analyze.

    Returns:
    pandas.Series: Median sea level change for the location relative to 2020.
    """
    sea_level_data = ds.sea_level_change.sel(locations=location_id)
    median_sea_level = sea_level_data.median(dim='samples')

    if use_2020_as_reference:
        sea_level_2020 = median_sea_level.values[0]
    else:
        sea_level_2020 = 0
    sea_level_change = (median_sea_level - sea_level_2020) / 1000

    return sea_level_change.to_series()

def process_location(datasets, location_id, scenario, output_dir, use_2020_as_reference=True, force_reprocess=False):
    """
    Process sea level change for a single location and save results.

    Args:
    datasets (dict): A dictionary of xarray.Datasets keyed by workflow ID.
    location_id (int): The ID of the location to analyze.
    scenario (str): The SSP scenario being analyzed.
    output_dir (str): Directory to save output CSV file.
    """
    output_file = os.path.join(output_dir, f'sea_level_change_loc{location_id}_{scenario}.csv')
    if not use_2020_as_reference:
        output_file = os.path.join(output_dir, f'sea_level_change_loc{location_id}_{scenario}_raw.csv')
    if os.path.exists(output_file) and not force_reprocess:
        print(f"Results for location {location_id} already saved to {output_file}, skipping processing.")
        return
    all_data = []
    for workflow_id, ds in datasets.items():
        print(f'Processing workflow_id: {workflow_id} for location_id: {location_id}')
        sea_level_change = get_median_sea_level_change(ds, location_id, use_2020_as_reference)
        sea_level_change = sea_level_change.reset_index()
        sea_level_change['workflow_id'] = workflow_id
        all_data.append(sea_level_change)

    result_df = pd.concat(all_data)

    result_df.to_csv(output_file, index=False)
    print(f"Results for location {location_id} saved to {output_file}")
    return result_df

def save_ar6_sea_level_projections(location_ids, scenario, output_dir, workflow_ids, use_2020_as_reference=True):
    """
    Main function to run the sea level projection analysis. The results are from the github repo.
    """
    datasets = load_ar6_data(scenario, workflow_ids)

    os.makedirs(output_dir, exist_ok=True)

    for location_id in location_ids:
        df = process_location(datasets, location_id, scenario, output_dir, use_2020_as_reference)
        print(df)

def rename_dir(dir, extension="csv"):
    """
    For all the csv files in the dir, add '_raw' to the filename
    """
    for file in os.listdir(dir):
        if file.endswith(f".{extension}"):
            os.rename(os.path.join(dir, file), os.path.join(dir, file.replace(f".{extension}", "_raw.csv")))

def save_data_from_ipcc_projection_tool_raw(ids: list|int, save_dir='data/ar6/projection_tool_raw'):
    """
    Save data from the IPCC projection tool.

    Args:
    ids (list): A list of PSMSL IDs for the locations.
    save_dir (str): The directory to save the downloaded Excel files.
    """
    if not isinstance(ids, list):
        ids = [ids]

    os.makedirs(save_dir, exist_ok=True)

    for id in ids:
        file_path = os.path.join(save_dir, f"ipcc_ar6_sea_level_projection_psmsl_id_{id}.xlsx")
        if os.path.exists(file_path):
            print(f"File {file_path} already exists, skipping download.")
            continue
        url = f"https://d3qt3aobtsas2h.cloudfront.net/edge/ws/search/projection?psmsl_id={id}&data_layer=scenario&format=csv"
        response = requests.get(url)

        if response.status_code == 200:
            with open(file_path, 'wb') as file:
                file.write(response.content)
            print(f"Data for ID {id} saved to {file_path}")
        else:
            print(f"Failed to download data for ID {id}. HTTP Status code: {response.status_code}")

def save_data_from_ipcc_projection_tool(ids: list|int,
                                        data_dir='data/ar6/projection_tool_raw',
                                        save_dir='data/ar6',
                                        **kwargs):
    """
    Save data from the IPCC projection tool, and process to csv format.
    The website uses confidence level to medium.

    Args:
    ids (list): A list of PSMSL IDs for the locations.
    data_dir (str): The directory to save the downloaded Excel files.
    save_dir (str): The directory to save the processed CSV files.
    kwargs:
        confidences (list or str): Confidence levels to filter.
        scenarios (list or str): Scenarios to filter.
        quantiles (list or int): Quantiles to filter.
    """
    save_data_from_ipcc_projection_tool_raw(ids, data_dir)
    os.makedirs(save_dir, exist_ok=True)
    files = list_files(data_dir)
    for file in files:
        df = process_ipcc_projection_tool_data(file, sheet_name='Total', **kwargs)
        df.to_csv(os.path.join(save_dir, os.path.basename(file).replace('.xlsx', '.csv')), index=False)
        print(f"Sea Level Projection Data for ID {id} saved to {os.path.join(save_dir, os.path.basename(file).replace('.xlsx', '.csv'))}")

        df = process_ipcc_projection_tool_data(file, sheet_name='VerticalLandMotion', **kwargs)
        save_path = os.path.join(save_dir, os.path.basename(file).replace('sea_level_projection', 'vlm_projection').replace('.xlsx', '.csv'))
        df.to_csv(save_path, index=False)
        print(f"VLM projection Data for ID {id} saved to {save_path}")

def list_files(dir, extension="xlsx"):
    """
    List all the files in the dir with the given extension
    """
    return [os.path.join(dir, file) for file in os.listdir(dir) if file.endswith(f".{extension}")]

def process_ipcc_projection_tool_data(excel_path, confidences=None, scenarios=None, quantiles=None, sheet_name='Total', interpolate_yearly: bool = True):
    """
    Process the IPCC projection data and return a formatted DataFrame.

    Parameters:
    - excel_path (str): Path to the Excel file.
    - confidences (list or str): Confidence levels to filter. Default is ['medium', 'low']
    - scenarios (list or str): Scenarios to filter. Default is ['ssp245', 'ssp126']
    - quantiles (list or int): Quantiles to filter. Default is [50]

    Returns:
    - pd.DataFrame: Formatted DataFrame with columns 'Year' and filtered data.
    """
    if confidences is None:
        confidences = ['medium', 'low']
    if scenarios is None:
        scenarios = ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]
    if quantiles is None:
        quantiles = [5, 50, 95]
    df = pd.read_excel(excel_path, sheet_name=sheet_name)

    if isinstance(confidences, str):
        confidences = [confidences]
    if isinstance(scenarios, str):
        scenarios = [scenarios]
    if isinstance(quantiles, int):
        quantiles = [quantiles]

    filtered_df = df[df['confidence'].isin(confidences) & df['scenario'].isin(scenarios) & df['quantile'].isin(quantiles)]

    years = list(range(2020, 2151, 10))
    result = {'Year': years}

    for index, row in filtered_df.iterrows():
        col_name = f"{row['scenario']}_{row['confidence']}_{row['quantile']}"
        result[col_name] = row[years].values

    result = pd.DataFrame(result)

    if interpolate_yearly:
        result = interpolate_yearly_data(result)
    return result



def interpolate_yearly_data(processed_data, year_col='Year'):
    """
    Interpolates the sea level change data to a yearly frequency.
    Not implemented yet.

    Parameters:
    - processed_data (pd.DataFrame): DataFrame containing the processed data with 10-year intervals.
    - year_col (str): Column name indicating the year. Default is 'Year'.

    Returns:
    - pd.DataFrame: DataFrame with interpolated yearly data.
    """
    years = processed_data[year_col]
    data_without_years = processed_data.drop(columns=[year_col])

    new_years = range(years.min(), years.max() + 1)
    interpolated_data = pd.DataFrame({year_col: new_years})

    for col in data_without_years.columns:
        f = interp1d(years, data_without_years[col], kind='linear', fill_value='extrapolate')
        interpolated_data[col] = f(new_years)

    return interpolated_data

def station2lonlat(station_name: str|int, data_dir="data/plots"):
    """
    Returns a dictionary with the longitude and latitude of the station.

    Parameters:
    station_name (str|int): The name of the station.
    data_dir (str): The directory of the data.

    Returns:
    dict: A dictionary with two elements: 'lon' and 'lat'.
    """
    possible_station_names = os.listdir(data_dir)
    if station_name in possible_station_names:
        station_info = pd.read_csv(os.path.join(data_dir, station_name, "GMR.csv"))
        return {"lon": station_info['lon'].iloc[0], "lat": station_info['lat'].iloc[0]}

    return get_ar6_station_latlon(station_name)

def ensure_ar6_location_list(file_path="data/ar6/location_list.txt"):
    """
    Ensure the AR6 location list file exists locally. If not, download it from the official GitHub repository.

    Args:
        file_path (str): Path to save the location list file.
    """
    url = "https://zenodo.org/records/6382554/files/location_list.lst?download=1"
    if not os.path.exists(file_path):
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        print(f"Downloading AR6 location list from {url} to {file_path} ...")
        response = requests.get(url)
        if response.status_code == 200:
            with open(file_path, "wb") as f:
                f.write(response.content)
            print("Download complete.")
        else:
            raise RuntimeError(f"Failed to download location list. Status code: {response.status_code}")
    else:
        print(f"Location list already exists at {file_path}.")

def read_ar6_location_list(file_path="data/ar6/location_list.txt"):
    """
    Read the location list from a file
    :param file_path: str, the path to the location list file
    :return: DataFrame with columns: station_name, station_id, lat, lon
    """
    ensure_ar6_location_list()
    return pd.read_csv(file_path, delimiter="\t", names=["station_name", "station_id", "lat", "lon"])

def get_ar6_station_id(station_name):
    location_list = read_ar6_location_list()
    station_info = location_list[location_list['station_name'] == station_name]
    if station_info.empty:
        raise ValueError(f"Station not found in the location list: {station_name}")
    return station_info['station_id'].values[0]

def get_ar6_station_name(station_id):
    location_list = read_ar6_location_list()
    station_info = location_list[location_list['station_id'] == station_id]
    if station_info.empty:
        raise ValueError(f"Station not found in the location list: {station_id}")
    return station_info['station_name'].values[0]

def get_ar6_station_info(station_id):
    """
    Get the latitude and longitude of a station
    :param station_id: str or int, the station id or the station name

    """
    location_list = read_ar6_location_list()

    if isinstance(station_id, str):
        station_info = location_list[location_list['station_name'] == station_id]
    else:
        station_info = location_list[location_list['station_id'] == station_id]

    if station_info.empty:
        raise ValueError(f"Station not found in the location list: {station_id}")

    return station_info

def get_ar6_station_latlon(station_id):
    """
    Get the latitude and longitude of a station
    :param station_id: str or int, the station id or the station name
    :return: dict with elements: lat, lon
    """
    station_info = get_ar6_station_info(station_id)
    return {'lat': station_info['lat'].values[0], 'lon': station_info['lon'].values[0]}

def get_station_name(station_id):
    """
    Get the name of a station
    :param station_id: int, the id of the station
    :return: str, the name of the station
    """
    if isinstance(station_id, (int, float)):
        station_id = int(station_id)
    else:
        raise ValueError(f"station_id must be an integer or float: {station_id}")

    for name, id in station_loc_id_map.items():
        if id == station_id:
            return name

    if station_id in station_loc_id_map:
        return station_id
    return get_ar6_station_info(station_id)['station_name'].values[0]

def get_station_id(station_name):
    """
    Get the id of a station
    :param station_name: str, the name of the station
    :return: int, the id of the station
    """
    try:
        return int(station_name)
    except ValueError:
        pass
    if isinstance(station_name, (int, float)):
        return int(station_name)

    if station_name in station_loc_id_map:
        return station_loc_id_map[station_name]
    return get_ar6_station_id(station_name)

def find_closest_ar6_station_latlon(lat, lon, location_list=None, max_distance=300, top_k=1):
    """
    Find the closest station based on latitude and longitude.

    Parameters:
    lat (float): The latitude of the point of interest.
    lon (float): The longitude of the point of interest.
    location_list (pd.DataFrame): The location list with columns: station_name, station_id, lat, lon.
    max_distance (float): Maximum allowed distance.

    Returns:
    int or list: The station_id(s) of the closest station(s).
    """
    if location_list is None:
        location_list = read_location_list()
    location_list['distance'] = np.sqrt((location_list['lat'] - lat)**2 + (location_list['lon'] - lon)**2)

    closest_stations = location_list.nsmallest(top_k, 'distance')

    if top_k == 1:
        if closest_stations.iloc[0]['distance'] > max_distance:
            raise ValueError("No station found within the max distance")
        return closest_stations.iloc[0]['station_id']
    else:
        valid_stations = closest_stations[closest_stations['distance'] <= max_distance]
        if valid_stations.empty:
            raise ValueError("No stations found within the max distance")
        return valid_stations['station_id'].tolist()

def find_closest_ar6_station(station, location_list=None, max_distance=2):
    """
    Find the closest AR6 station based on the station name.

    Parameters:
    station (str): The name of the station. Only old station names are supported.
    location_list (pd.DataFrame): The location list with columns: station_name, station_id, lat, lon.
    max_distance (float): Maximum allowed distance.

    Returns:
    int: The station_id of the closest station.
    """
    if location_list is None:
        location_list = read_location_list()

    loc_id = station_loc_id_map.get(station)
    if loc_id is not None:
        return loc_id

    station_coords = station2lonlat(station)
    return find_closest_ar6_station_latlon(station_coords['lat'], station_coords['lon'], location_list, max_distance)


def read_location_list(file_path="data/ar6/location_list.txt"):
    """
    Read the location list from a file.

    Parameters:
    file_path (str): The path to the location list file.

    Returns:
    pd.DataFrame: DataFrame with columns: station_name, station_id, lat, lon.
    """
    return pd.read_csv(file_path, delimiter="\t", names=["station_name", "station_id", "lat", "lon"])

if __name__ == "__main__":
    df = read_ar6_location_list()
    location_ids = df[df['station_id'] < 2358]['station_id'].tolist()
    save_data_from_ipcc_projection_tool(location_ids)








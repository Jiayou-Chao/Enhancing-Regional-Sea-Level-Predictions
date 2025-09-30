import cdsapi
import os

c = cdsapi.Client()

base_dir = "data/era5data/"
vars = ["10m_u_component_of_wind", "10m_v_component_of_wind", "2m_dewpoint_temperature",
"2m_temperature", "surface_pressure", "total_precipitation"]
for var in vars:
    os.makedirs(base_dir + var, exist_ok=True)

    for y in range(1950, 2025):
        year = str(y)
        filename = os.path.join(base_dir, var, f"{var}_{year}.grib")
        if os.path.exists(filename):
            print(f"File {filename} already exists, skipping download.")
            continue

        c.retrieve(
            'reanalysis-era5-land-monthly-means',
            {
                'product_type': 'monthly_averaged_reanalysis',
                'variable': var,
                'year': [
                    year
                ],
                'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                ],
                'time': '00:00',
                'format': 'grib',
            },
            filename)

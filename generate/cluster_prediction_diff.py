import pandas as pd
import hdbscan
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import argparse
from pathlib import Path
import json
import numpy as np

def load_and_preprocess_data(file_path: str) -> pd.DataFrame:
    """
    Loads station data from a CSV, handles missing values, and scales
    the relevant features for clustering.

    Args:
        file_path: The path to the input CSV file.
                   The file should contain 'lat', 'lon', 'scenario1', and 'scenario2' columns.

    Returns:
        A DataFrame with the original data and scaled columns for clustering.
        Returns an empty DataFrame if the file is not found or is empty.
    """
    if not Path(file_path).is_file():
        print(f"Error: File not found at {file_path}")
        return pd.DataFrame()

    data = pd.read_csv(file_path)

    data.dropna(subset=['lat', 'lon', 'scenario1', 'scenario2'], inplace=True)

    if data.empty:
        print("Warning: Data is empty after dropping NA values.")
        return pd.DataFrame()

    features_to_scale = ['lat', 'lon', 'diff']
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data[features_to_scale])

    data['scaled_lat'] = data_scaled[:, 0]
    data['scaled_lon'] = data_scaled[:, 1]
    data['scaled_diff'] = data_scaled[:, 2]

    return data

def perform_hdbscan_clustering(data: pd.DataFrame) -> pd.DataFrame:
    """
    Performs HDBSCAN clustering on the preprocessed data.

    Args:
        data: A DataFrame containing the scaled 'lat', 'lon', and 'diff'
              columns, typically from `load_and_preprocess_data`.

    Returns:
        The input DataFrame with an added 'cluster_id' column.
    """
    features = data[['scaled_lat', 'scaled_lon', 'scaled_diff']]

    clusterer = hdbscan.HDBSCAN(min_cluster_size=5, gen_min_span_tree=True)
    clusterer.fit(features)

    data['cluster_id'] = clusterer.labels_

    return data

def visualize_clusters(data: pd.DataFrame, output_path: str):
    """
    Generates and saves an interactive map of the clusters using Plotly.

    Args:
        data: The DataFrame containing station data with a 'cluster_id' column.
        output_path: The path to save the HTML file for the plot.
    """
    plot_data = data.copy()
    unique_ids = sorted(plot_data['cluster_id'].unique())
    category_order = [str(uid) for uid in unique_ids]

    plot_data['cluster_id'] = plot_data['cluster_id'].astype(str)

    fig = px.scatter_geo(
        plot_data,
        lat='lat',
        lon='lon',
        color='cluster_id',
        hover_name='station_id',
        hover_data={'lat': True, 'lon': True, 'diff': True, 'cluster_id': True},
        title='HDBSCAN Clustering of Stations',
        projection='natural earth',
        category_orders={"cluster_id": category_order}
    )

    fig.write_html(output_path)
    print(f"Cluster visualization saved to {output_path}")


def save_clusters_static_image(data: pd.DataFrame, output_path: str):
    """
    Saves a static PNG image of the cluster visualization for use in publications.

    Args:
        data: DataFrame with station data and a 'cluster_id' column.
        output_path: Path to save the PNG image.
    """
    plot_data = data.copy()
    unique_ids = sorted(plot_data['cluster_id'].unique())
    category_order = [str(uid) for uid in unique_ids]

    plot_data['cluster_id'] = plot_data['cluster_id'].astype(str)

    fig = px.scatter_geo(
        plot_data,
        lat='lat',
        lon='lon',
        color='cluster_id',
        hover_name='station_id',
        hover_data={'lat': True, 'lon': True, 'diff': True, 'cluster_id': True},
        title='HDBSCAN Clustering of Stations',
        projection='natural earth',
        category_orders={"cluster_id": category_order}
    )
    fig.update_layout(showlegend=True)
    fig.write_image(output_path, format="png", scale=3)
    print(f"Static cluster visualization saved to {output_path}")

def get_descriptive_statistics(series: pd.Series) -> dict:
    """
    Computes a standardized set of descriptive statistics for a pandas Series.
    """
    return {
        'mean': series.mean(),
        'std': series.std(),
        'min': series.min(),
        'max': series.max(),
        'p25': series.quantile(0.25),
        'p50': series.quantile(0.50),
        'p75': series.quantile(0.75),
    }

def analyze_and_save_statistics(data: pd.DataFrame, output_path: str):
    """
    Calculates comprehensive statistics for each cluster and saves them to a JSON file.
    """
    cluster_ids = sorted(data['cluster_id'].unique())

    data['abs_diff'] = (data['scenario1'] - data['scenario2']).abs()

    with np.errstate(divide='ignore', invalid='ignore'):
        percentage_increase = (data['scenario1'] - data['scenario2']) / data['scenario2']

    data['percentage_increase'] = percentage_increase.replace([np.inf, -np.inf], None)

    overall_stats = {
        "count": len(data),
        "location": {
            "mean_lat": data['lat'].mean(),
            "mean_lon": data['lon'].mean()
        },
        "diff": get_descriptive_statistics(data['diff']),
        "abs_diff": get_descriptive_statistics(data['abs_diff']),
        "percentage_increase": get_descriptive_statistics(data['percentage_increase'].dropna()),
        "mse": (data['diff']**2).mean()
    }

    cluster_stats = {}
    for cluster_id in cluster_ids:
        cluster_data = data[data['cluster_id'] == cluster_id].copy()

        with np.errstate(divide='ignore', invalid='ignore'):
            cluster_percentage_increase = (cluster_data['scenario1'] - cluster_data['scenario2']) / cluster_data['scenario2']
        cluster_data['percentage_increase'] = cluster_percentage_increase.replace([np.inf, -np.inf], None)

        cluster_stats[str(cluster_id)] = {
            "count": len(cluster_data),
            "location": {
                "mean_lat": cluster_data['lat'].mean(),
                "mean_lon": cluster_data['lon'].mean()
            },
            "diff": get_descriptive_statistics(cluster_data['diff']),
            "abs_diff": get_descriptive_statistics(cluster_data['abs_diff']),
            "percentage_increase": get_descriptive_statistics(cluster_data['percentage_increase'].dropna()),
            "mse": (cluster_data['diff']**2).mean()
        }

    final_json = {
        "overall": overall_stats,
        "clusters": cluster_stats
    }

    with open(output_path, 'w') as f:
        json.dump(final_json, f, indent=4, allow_nan=False)

    print(f"Cluster statistics saved to {output_path}")

def main():
    """
    Main execution function to run the clustering and analysis pipeline.
    """
    parser = argparse.ArgumentParser(
        description="Perform HDBSCAN clustering on tidal gauge station data."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input CSV file containing station data."
    )
    parser.add_argument(
        "-o", "--output_dir",
        type=str,
        default="results/clusters",
        help="Directory to save output files (default: results/clusters)"
    )
    parser.add_argument(
        "--min_cluster_size",
        type=int,
        default=5,
        help="Minimum cluster size for HDBSCAN (default: 5)"
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    visualization_path = output_dir / "cluster_visualization.html"
    statistics_path = output_dir / "cluster_statistics.json"
    clustered_csv_path = output_dir / "clustered_data.csv"
    static_image_path = output_dir / "cluster_visualization.png"

    print(f"Loading data from {args.input_file}...")
    data = load_and_preprocess_data(args.input_file)

    if data.empty:
        print("Error: No valid data to process. Exiting.")
        return

    print("Performing HDBSCAN clustering...")
    clustered_data = perform_hdbscan_clustering(data)

    clustered_data.to_csv(clustered_csv_path, index=False)
    print(f"Clustered data saved to {clustered_csv_path}")

    print("Generating visualization...")
    visualize_clusters(clustered_data, str(visualization_path))

    print("Saving static image of the clusters...")
    save_clusters_static_image(clustered_data, str(static_image_path))

    print("Analyzing cluster statistics...")
    analyze_and_save_statistics(clustered_data, str(statistics_path))

    print("\nClustering complete!")
    print(f"- Clustered data saved to: {clustered_csv_path}")
    print(f"- Interactive visualization saved to: {visualization_path}")
    print(f"- Static image of clusters saved to: {static_image_path}")
    print(f"- Statistics saved to: {statistics_path}")

if __name__ == "__main__":
    main()
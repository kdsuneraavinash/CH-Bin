from configparser import SectionProxy
from pathlib import Path

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from bin_x.cli.utils import handle_error, reduce_dimensions_to_2d
from bin_x.core.clustering.algorithm import fit_cluster
from bin_x.core.clustering.distance_matrix import create_distance_matrix
from bin_x.core.config import USER_CONFIG


def _visualize_final_result(
    num_clusters: int, df_features: pd.DataFrame, df_dist_bin: pd.DataFrame, split_bin_png: Path, bin_counts_png: Path
):
    df_kmer_bins = df_dist_bin.merge(df_features, left_on="CONTIG_NAME", right_on="PARENT_NAME")
    df_kmer_counts = df_kmer_bins.drop(["CONTIG_NAME_x", "CONTIG_NAME_y", "PARENT_NAME", "CLUSTER", "BIN"], axis=1)
    df_2d = reduce_dimensions_to_2d(df_kmer_counts)

    plt.figure(figsize=(15, 10))
    plt.title("Convex hull Result (Split)")
    for i in range(-1, num_clusters):
        sns.scatterplot(x="X", y="Y", data=df_2d[df_kmer_bins.BIN == i], label=f"Cluster {i}", alpha=0.1)
    plt.savefig(split_bin_png)
    plt.close()

    plt.figure(figsize=(15, 10))
    plt.title("Bin Counts")
    sns.histplot(df_kmer_bins, x="BIN", discrete=True)
    plt.savefig(bin_counts_png)
    plt.close()


def perform_clustering(
    features_csv: Path,
    operating_dir: Path,
    num_neighbors: int = 15,
    max_iterations: int = 10,
    metric: str = "convex",
    qp_solver: str = "quadprog",
) -> Path:
    """
    Perform binning and output the binning result.

    :param features_csv: CSV containing the feature vectors and initial bins.
    :param operating_dir: Directory to write temp files to.
    :param num_neighbors: Number of neighbors to consider for polytope.
    :param max_iterations: Number of maximum iterations to perform.
    :param metric: Polytope distance matrix (convex/affine)
    :param qp_solver: Quadratic programming problem solver. (quadprog/cvxopt)
    :return: Path of the binning result dataset.
    """
    dist_bin_csv = operating_dir / "bin.csv"
    split_bin_png = operating_dir / "split_bin.png"
    bin_counts_png = operating_dir / "bin_counts.png"
    operating_dir.mkdir(parents=True, exist_ok=True)

    # 01. Read feature CSV
    click.secho("01. Reading feature CSV...", bold=True)
    df_features = pd.read_csv(features_csv)

    num_clusters = df_features.CLUSTER.max() + 1
    initial_bins: np.ndarray = df_features.CLUSTER.values.copy()
    samples: np.ndarray = df_features.drop(["CONTIG_NAME", "PARENT_NAME", "CLUSTER"], axis=1).values
    num_samples = len(samples)

    # 02. Create a distance matrix
    click.secho(f"02. Creating a distance matrix of {num_samples}x{num_samples} shape...", bold=True)
    distance_matrix = create_distance_matrix(samples, operating_dir)

    # 03. Perform binning using specified solver and metric
    click.secho(f"03. Performing binning using {qp_solver} solver", bold=True)
    convex_labels = fit_cluster(
        num_samples=num_samples,
        samples=samples,
        num_clusters=num_clusters,
        distance_matrix=distance_matrix,
        initial_bins=initial_bins,
        num_neighbors=num_neighbors,
        max_iterations=max_iterations,
        metric=metric,
        qp_solver=qp_solver,
    )
    if np.any(convex_labels < 0):
        raise ValueError("There were some un-clustered points left... Aborting.")

    # 04. Assigning bins with majority voting (If there were more than one voting column)
    click.secho("04. Assigning bins", bold=True)
    df_samples: pd.DataFrame = df_features.drop("CLUSTER", axis=1)
    df_bin_column: pd.DataFrame = pd.DataFrame({"BIN": convex_labels})

    df_combined: pd.DataFrame = pd.concat([df_samples, df_bin_column], axis=1)
    parent_groups = df_combined[["PARENT_NAME", "BIN"]].groupby("PARENT_NAME")
    df_dist_bin: pd.DataFrame = parent_groups.BIN.apply(lambda x: np.bincount(x).argmax()).reset_index()
    df_dist_bin.rename(columns={"PARENT_NAME": "CONTIG_NAME"}, inplace=True)
    df_dist_bin.to_csv(dist_bin_csv, index=False)  # noqa
    click.secho(f"Dumped binning assignment CSV at {dist_bin_csv}", fg="green", bold=True)

    # 05. Visualize
    click.secho("08. Drawing plots...", bold=True)
    _visualize_final_result(num_clusters, df_features, df_dist_bin, split_bin_png, bin_counts_png)

    return dist_bin_csv


def run_perform_clustering(features_csv: Path, operating_dir: Path, parameters: SectionProxy) -> Path:
    """
    Perform binning and output the binning result.

    :param features_csv: CSV containing the feature vectors and initial bins.
    :param operating_dir: Directory to write temp files to.
    :param parameters: Parameters INI section.
    :return: Path of the binning result dataset.
    """

    return perform_clustering(
        features_csv=features_csv,
        operating_dir=operating_dir,
        num_neighbors=int(parameters["AlgoNumNeighbors"]),
        max_iterations=int(parameters["AlgoMaxIterations"]),
        metric=parameters["AlgoDistanceMetric"],
        qp_solver=parameters["AlgoQpSolver"],
    )


@click.command()
@click.option("--config", prompt="Configuration file", help="The INI File to use for tool configuration.", type=Path)
@click.option("--features", prompt="Features file", help="The features CSV file with initial clustering.", type=Path)
@click.option("--out", prompt="Output Directory", help="The output directory for the tool.", type=Path)
def main(config: Path, features: Path, out: Path):
    try:
        USER_CONFIG.read(config)
        parameters = USER_CONFIG["PARAMETERS"]
        run_perform_clustering(features, out, parameters)
    except Exception as e:
        handle_error(e)


if __name__ == "__main__":
    main()

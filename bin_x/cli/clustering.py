import os
from configparser import SectionProxy
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
from numpy.lib.format import open_memmap  # noqa

from bin_x.core.clustering.algorithm import fit_cluster
from bin_x.core.clustering.distance_matrix import (
    create_distance_matrix,
    create_in_mem_distance_matrix,
)
from bin_x.core.clustering.dump_bins import dump_bins


def perform_clustering(
    contig_fasta: Path,
    features_csv: Path,
    operating_dir: Path,
    num_neighbors: int = 15,
    max_iterations: int = 10,
    metric: str = "convex",
    qp_solver: str = "quadprog",
    in_mem_dist_matrix: bool = True,
    distance_matrix_cache: Optional[Path] = None,
) -> Path:
    """
    Perform binning and output the binning result.

    :param contig_fasta: Contig file used for feature generation.
    :param features_csv: CSV containing the feature vectors and initial bins.
    :param operating_dir: Directory to write temp files to.
    :param num_neighbors: Number of neighbors to consider for polytope.
    :param max_iterations: Number of maximum iterations to perform.
    :param metric: Polytope distance matrix (convex/affine)
    :param qp_solver: Quadratic programming problem solver. (quadprog/cvxopt)
    :param in_mem_dist_matrix: Whether to use in memory distance matrix.
    :param distance_matrix_cache: Previously computed distance matrix.
                                    in_mem_dist_matrix must be false if this is provided.
    :return: Path of the binning result dataset.
    """
    dist_bin_csv = operating_dir / "binning-assignment.csv"
    bin_dump_dir = operating_dir / "bins"
    operating_dir.mkdir(parents=True, exist_ok=True)
    bin_dump_dir.mkdir(parents=True, exist_ok=True)

    # 01. Read feature CSV
    click.secho(">> Reading feature CSV...", fg="green", bold=True)
    df_features = pd.read_csv(features_csv)

    num_clusters = df_features.CLUSTER.max() + 1
    initial_bins: np.ndarray = df_features.CLUSTER.values.copy()
    samples: np.ndarray = df_features.drop(["CONTIG_NAME", "PARENT_NAME", "CLUSTER"], axis=1).values
    num_samples = len(samples)

    # 02. Create a distance matrix
    distance_matrix_filename = distance_matrix_cache
    if in_mem_dist_matrix:
        assert distance_matrix_filename is None, "Distance matrix cache cannot be used in in-memory mode"
        click.secho(f">> Creating the distance matrix({num_samples}x{num_samples}) in-memory...", fg="green", bold=True)
        distance_matrix = create_in_mem_distance_matrix(samples)
    else:
        if distance_matrix_filename is None:
            click.secho(f">> Creating a distance matrix of {num_samples}x{num_samples} shape...", fg="green", bold=True)
            distance_matrix_filename = create_distance_matrix(samples, operating_dir)
        else:
            click.secho(f">> Reusing distance matrix at {distance_matrix_filename}...", fg="green", bold=True)
        distance_matrix = open_memmap(filename=distance_matrix_filename, mode="r", shape=(num_samples, num_samples))

    # 03. Perform binning using specified solver and metric
    click.secho(f">> Performing binning using {qp_solver} solver", fg="green", bold=True)
    convex_labels = fit_cluster(
        samples=samples,
        num_clusters=num_clusters,
        distance_matrix=distance_matrix,
        initial_bins=initial_bins,
        num_neighbors=num_neighbors,
        max_iterations=max_iterations,
        metric=metric,
        qp_solver=qp_solver,
    )

    # Delete the distance matrix file (dont delete if using a cache)
    if not in_mem_dist_matrix and distance_matrix_cache is None:
        os.remove(distance_matrix_filename)
    if np.any(convex_labels < 0):
        raise ValueError("There were some un-clustered points left... Aborting.")

    # 04. Assigning bins with majority voting (If there were more than one voting column)
    click.secho(">> Assigning bins", fg="green", bold=True)
    df_samples: pd.DataFrame = df_features.drop("CLUSTER", axis=1)
    df_bin_column: pd.DataFrame = pd.DataFrame({"BIN": convex_labels})

    df_combined: pd.DataFrame = pd.concat([df_samples, df_bin_column], axis=1)
    parent_groups = df_combined[["PARENT_NAME", "BIN"]].groupby("PARENT_NAME")
    df_dist_bin: pd.DataFrame = parent_groups.BIN.apply(lambda x: np.bincount(x).argmax()).reset_index()
    df_dist_bin.rename(columns={"PARENT_NAME": "CONTIG_NAME"}, inplace=True)
    df_dist_bin.to_csv(dist_bin_csv, index=False)  # noqa
    click.secho(f"Dumped binning assignment CSV at {dist_bin_csv}", bold=True)

    # 05. Creating bin files with contigs
    click.secho(">> Writing binned FASTA files", fg="green", bold=True)
    dump_bins(df_dist_bin, contig_fasta, bin_dump_dir)
    click.secho(f"Dumped binned fasta to {bin_dump_dir}", bold=True)

    return dist_bin_csv


def run_perform_clustering(
    contig_fasta: Path,
    features_csv: Path,
    operating_dir: Path,
    parameters: SectionProxy,
    distance_matrix_cache: Optional[Path] = None,
) -> Path:
    """
    Perform binning and output the binning result.

    :param contig_fasta: Contig file used for feature generation.
    :param features_csv: CSV containing the feature vectors and initial bins.
    :param operating_dir: Directory to write temp files to.
    :param parameters: Parameters INI section.
    :param distance_matrix_cache: Previously computed distance matrix.
    :return: Path of the binning result dataset.
    """

    return perform_clustering(
        contig_fasta=contig_fasta,
        features_csv=features_csv,
        operating_dir=operating_dir,
        num_neighbors=int(parameters["AlgoNumNeighbors"]),
        max_iterations=int(parameters["AlgoMaxIterations"]),
        metric=parameters["AlgoDistanceMetric"],
        qp_solver=parameters["AlgoQpSolver"],
        in_mem_dist_matrix=parameters.getboolean("InMemDistMatrix"),
        distance_matrix_cache=distance_matrix_cache,
    )

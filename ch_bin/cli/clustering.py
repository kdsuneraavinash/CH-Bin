import logging
from configparser import SectionProxy
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.lib.format import open_memmap  # noqa

from ch_bin.core.clustering.algorithm import fit_cluster
from ch_bin.core.clustering.distance_matrix import (
    create_distance_matrix,
    create_in_mem_distance_matrix,
)
from ch_bin.core.clustering.dump_bins import dump_bins

logger = logging.getLogger(__name__)


def perform_clustering(
    contig_fasta: Path,
    features_csv: Path,
    operating_dir: Path,
    num_neighbors: int = 15,
    max_iterations: int = 10,
    metric: str = "convex",
    qp_solver: str = "quadprog",
    in_mem_dist_matrix: bool = True,
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
    :return: Path of the binning result dataset.
    """
    dist_bin_csv = operating_dir / "binning-assignment.csv"
    bin_dump_dir = operating_dir / "bins"
    operating_dir.mkdir(parents=True, exist_ok=True)
    bin_dump_dir.mkdir(parents=True, exist_ok=True)

    # 01. Read feature CSV
    logger.info(">> Reading feature CSV...")
    df_features = pd.read_csv(features_csv)

    num_clusters = df_features.CLUSTER.max() + 1
    initial_bins: np.ndarray = df_features.CLUSTER.values.copy()
    samples: np.ndarray = df_features.drop(["CONTIG_NAME", "PARENT_NAME", "CLUSTER"], axis=1).values
    num_samples = len(samples)

    # 02. Create a distance matrix
    if in_mem_dist_matrix:
        logger.info(">> Creating a distance matrix of %s in-memory...", (num_samples, num_samples))
        distance_matrix = create_in_mem_distance_matrix(samples)
    else:
        logger.info(">> Creating a distance matrix of %s in-disk...", (num_samples, num_samples))
        distance_matrix_filename = create_distance_matrix(samples, operating_dir)
        distance_matrix = open_memmap(filename=distance_matrix_filename, mode="r", shape=(num_samples, num_samples))

    # 03. Perform binning using specified solver and metric
    logger.info(">> Performing binning using %s solver...", qp_solver)
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
    if np.any(convex_labels < 0):
        raise ValueError("There were some un-clustered points left... Aborting.")

    # 04. Assigning bins with majority voting (If there were more than one voting column)
    logger.info(">> Assigning bins...")
    df_samples: pd.DataFrame = df_features.drop("CLUSTER", axis=1)
    df_bin_column: pd.DataFrame = pd.DataFrame({"BIN": convex_labels})

    df_combined: pd.DataFrame = pd.concat([df_samples, df_bin_column], axis=1)
    parent_groups = df_combined[["PARENT_NAME", "BIN"]].groupby("PARENT_NAME")
    df_dist_bin: pd.DataFrame = parent_groups.BIN.apply(lambda x: np.bincount(x).argmax()).reset_index()
    df_dist_bin.rename(columns={"PARENT_NAME": "CONTIG_NAME"}, inplace=True)
    df_dist_bin.to_csv(dist_bin_csv, index=False)  # noqa
    logger.info("Dumped binning assignment CSV at %s...", dist_bin_csv)

    # 05. Creating bin files with contigs
    logger.info(">> Writing binned FASTA files...")
    dump_bins(df_dist_bin, contig_fasta, bin_dump_dir)
    logger.info("Dumped binned FASTA files to %s...", bin_dump_dir)

    return dist_bin_csv


def run_perform_clustering(
    contig_fasta: Path,
    features_csv: Path,
    operating_dir: Path,
    parameters: SectionProxy,
) -> Path:
    """
    Perform binning and output the binning result.

    :param contig_fasta: Contig file used for feature generation.
    :param features_csv: CSV containing the feature vectors and initial bins.
    :param operating_dir: Directory to write temp files to.
    :param parameters: Parameters INI section.
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
    )

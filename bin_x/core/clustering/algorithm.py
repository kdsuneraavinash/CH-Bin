import numpy as np
from tqdm import tqdm

from bin_x.core.clustering.distance_matrix import find_m_nearest_neighbors
from bin_x.core.clustering.hull_distance import (
    affine_hull_distance,
    convex_hull_distance,
)


def _calculate_distance(x: np.ndarray, mat_p: np.ndarray, qp_solver: str, metric: str) -> float:
    """
    Function to find the distance from point x to a polytope defined by set P.
    Uses the quadratic optimization to find the distance.

    :param x: Query point.
    :param mat_p: Matrix with polytope points.
    :param qp_solver: Quadratic program solver algorithm.
    :param metric: Metric to use.
    :return: The distance from the query point to the polytope.
    """

    if metric == "convex":
        return convex_hull_distance(x, mat_p, solver=qp_solver)
    if metric == "affine":
        return affine_hull_distance(x, mat_p, solver=qp_solver)
    raise NotImplementedError(f"Metric {metric} not implemented")


def fit_cluster(
    num_samples: int,
    samples: np.ndarray,
    num_clusters: int,
    initial_bins: np.ndarray,
    distance_matrix: np.ndarray,
    num_neighbors: int = 15,
    max_iterations: int = 10,
    metric: str = "convex",
    qp_solver: str = "quadprog",
) -> np.ndarray:
    """
    Perform the Cevikalp et. al. 2019 convex hull binning algorithm
    specialized for metagenomic binning.

    :param num_samples: Number of samples in given dataset.
    :param samples: Dataset points.
    :param num_clusters: Number of clusters.
    :param num_neighbors: Number of neighbors to consider for polytope.
    :param distance_matrix: Pre-calculated distance matrix.
    :param initial_bins: Initial bin vector. Use -1 for un-binned.
    :param max_iterations: Number of maximum iterations to perform.
    :param metric: Polytope distance matrix (convex/affine)
    :param qp_solver: Quadratic programming problem solver. (quadprog/cvxopt)
    :return: Final binning result.
    """

    curr_bins: np.ndarray = initial_bins.copy()

    for i_iter in range(max_iterations):

        for i, i_sample in tqdm(
            enumerate(np.random.permutation(num_samples)), desc=f"Iteration {i_iter + 1}", total=num_samples, ncols=80
        ):
            min_distance: float = np.inf
            min_cluster: int = curr_bins[i_sample]

            # Reassign point to the closest cluster.
            curr_bins[i_sample] = -1
            neighbor_idx = find_m_nearest_neighbors(distance_matrix[i_sample], num_clusters, curr_bins, m=num_neighbors)
            if len(neighbor_idx) == 0:
                curr_bins[i_sample] = min_cluster

            for c in range(num_clusters):
                # Find convex hull distance
                distance = _calculate_distance(samples[i_sample], samples[neighbor_idx[c]], qp_solver, metric)
                if min_distance > distance:
                    min_distance, min_cluster = distance, c

            curr_bins[i_sample] = min_cluster  # Assign xi to cluster with smallest distance

        # If the assignments did not change, break
        diff_bins = initial_bins != curr_bins
        if not np.any(diff_bins):
            print("No changes done with previous iteration... Stopping at iteration", i_iter + 1)
            break
        print(f"Points changing cluster : {np.average(diff_bins) * 100}% ({sum(diff_bins)} points)")

        initial_bins = curr_bins
        curr_bins = curr_bins.copy()

    else:
        print("Exit due to max iteration limit")
    return curr_bins

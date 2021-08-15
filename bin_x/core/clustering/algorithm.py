import numpy as np
from tqdm import tqdm

from bin_x.core.clustering.distance_matrix import find_nearest_from_cluster
from bin_x.core.clustering.hull_distance import calculate_distance


def fit_cluster(
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
    points_to_assign: np.ndarray = np.where(curr_bins == -1)[0]
    num_points_to_assign = len(points_to_assign)

    distance_row = np.empty(shape=(len(samples),))
    for i_iter in range(max_iterations):

        sample_perm = np.random.permutation(points_to_assign)
        for i_sample in tqdm(sample_perm, desc=f"Iteration {i_iter + 1}", total=num_points_to_assign, ncols=80):
            min_distance: float = np.inf
            min_cluster: int = curr_bins[i_sample]

            curr_bins[i_sample] = -1  # Remove the point from current cluster
            distance_row[:] = distance_matrix[i_sample, :]
            for c in range(num_clusters):
                # 01. Find m closest points from the cluster
                # 02. Find convex hull distance
                cluster_point_idx = find_nearest_from_cluster(c, curr_bins, distance_row, num_neighbors)
                distance = calculate_distance(samples[i_sample], samples[cluster_point_idx], qp_solver, metric)
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

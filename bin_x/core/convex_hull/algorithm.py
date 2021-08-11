import time

import numpy as np

from bin_x.core.convex_hull.distance import distance_convex, find_m_nearest_idx


def _calculate_distance(x: np.ndarray, mat_p: np.ndarray, qp_solver: str, metric: str):
    """
    Function to find the distance from point x to a polytope defined by set P.
    Uses the quadratic optimization to find the distance.
    """

    if metric == "convex":
        return distance_convex(x, mat_p, solver=qp_solver)
    raise NotImplementedError(f"Metric {metric} not implemented")


def fit_cluster(
    n_iter: int,
    n_samples: int,
    samples: np.ndarray,
    n_clusters: int,
    n_neighbors: int,
    mat_dist: np.ndarray,
    init_bins: np.ndarray,
    metric: str = "convex",
    qp_solver: str = "quadprog",
) -> np.ndarray:
    curr_bins = init_bins.copy()

    for i_iter in range(n_iter):
        start_timestamp = time.time()

        for i, i_sample in enumerate(np.random.permutation(n_samples)):
            print(f"Iteration {i_iter + 1} progress: {i}/{n_samples}", end="\r")
            min_distance, min_cluster = np.inf, curr_bins[i_sample]

            for c in range(n_clusters):
                # Determine the m nearest samples of xi from cluster c
                (selected_idx,) = np.where(curr_bins == c)
                selected_idx = find_m_nearest_idx(n_neighbors, i_sample, mat_dist, selected_idx)
                # Find convex hull distance
                distance, _, __ = _calculate_distance(samples[i_sample], samples[selected_idx], qp_solver, metric)
                if min_distance > distance:
                    min_distance, min_cluster = distance, c

            curr_bins[i_sample] = min_cluster  # Assign xi to cluster with smallest distance

        end_timestamp = time.time()
        print(f"Iteration {i_iter + 1} took {end_timestamp - start_timestamp}s")

        # If the assignments did not change, break
        diff_bins = init_bins != curr_bins
        if not np.any(diff_bins):
            print("No changes done with previous iteration... " "Stopping at iteration", i_iter + 1)
            break
        print(f"Points changing cluster : {np.average(diff_bins) * 100}%")
        print()

        init_bins = curr_bins
        curr_bins = curr_bins.copy()

    else:
        print("Exit due to max iteration limit")
    return curr_bins

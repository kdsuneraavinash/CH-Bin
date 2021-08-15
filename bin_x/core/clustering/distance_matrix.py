import time
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap
from scipy.spatial.distance import cdist


def create_distance_matrix(arr: np.ndarray, operating_dir: Path) -> Path:
    """
    This calculates a distance matrix D of nxn size for given array of n elements.
    The (i,j) element of matrix will have the euclidean distance from ith point to the jth point.
    """
    n = len(arr)
    start_time = time.time()
    filename = operating_dir / "distance_matrix.npy"
    result = open_memmap(filename=filename, mode="w+", shape=(n, n))
    cdist(arr, arr, metric="euclidean", out=result)
    result.flush()
    print(f"Distance matrix calculated in {time.time() - start_time}s")
    return filename


def find_m_nearest_neighbors(distance_row: np.ndarray, n_clusters: int, clusters: np.ndarray, m: int):
    """
    Finds the indices of m nearest neighbors of the distance matrix row given.
    This will return all nearest neighbors in all clusters.
    The result will be a list of indices with indices being of the original array that created distance matrix.

    Note: p will not be considered as a nearest neighbor to p.

    :param distance_row: Number of nearest neighbors to find.
    :param n_clusters: Index of point to find the neighbors.
    :param clusters: Distance matrix.
    :param m: Filter for the required points.
    :return: M nearest neighbor indices.
    """
    result = []
    for c in range(n_clusters):
        idx = np.where(clusters == c)[0]
        if len(idx) <= m:
            result.append(idx)
        else:
            res_idx = np.argpartition(distance_row[idx], m)[:m]
            result.append(idx[res_idx])
    return result

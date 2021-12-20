import logging
import time
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)


def create_distance_matrix(arr: np.ndarray, operating_dir: Path) -> Path:
    """
    This calculates a distance matrix D of nxn size for given array of n elements.
    The (i,j) element of matrix will have the euclidean distance from ith point to the jth point.
    """
    n = len(arr)
    filename = operating_dir / "distance_matrix.npy"
    if filename.exists():
        logger.info("Reusing already existing distance matrix at %s.", filename)
        logger.debug("Assuming memmap shape %s", (n, n))
        return filename
    start_time = time.time()
    logger.debug("Started creating distance matrix at %s.", filename)
    result = open_memmap(filename=filename, mode="w+", shape=(n, n))
    cdist(arr, arr, metric="euclidean", out=result)
    result.flush()
    logger.debug("Ended creating distance matrix. Shape is %s", result.shape)
    logger.debug("Distance matrix calculated in %s s.", time.time() - start_time)
    return filename


def create_in_mem_distance_matrix(arr: np.ndarray) -> np.ndarray:
    """
    This calculates a distance matrix D of nxn size for given array of n elements.
    The (i,j) element of matrix will have the euclidean distance from ith point to the jth point.
    The created distance matrix would be in-memory.
    """
    start_time = time.time()
    logger.debug("Started creating distance matrix in-memory.")
    result = cdist(arr, arr, metric="euclidean")
    logger.debug("Ended creating distance matrix. Shape is %s", result.shape)
    logger.debug("Distance matrix calculated in %s s.", time.time() - start_time)
    return result


def find_nearest_from_cluster(c: int, curr_bins: np.ndarray, distance_row: np.ndarray, m: int) -> np.ndarray:
    """
    Find m closest points from the a cluster to a query sample point (denoted by its distance matrix row).

    :param c: Processing cluster number.
    :param curr_bins: Array with cluster assignments of all points.
    :param distance_row: Distance matrix row for the sample.
    :param m: Number of nearest points to consider.
    :return: All indices of nearest sample points.
    """
    cluster_point_idx = np.where(curr_bins == c)[0]
    if len(cluster_point_idx) <= m:
        return cluster_point_idx
    distance_values = distance_row[cluster_point_idx]
    closest_points_filter = np.argpartition(distance_values, kth=m - 1)[:m]
    return cluster_point_idx[closest_points_filter]

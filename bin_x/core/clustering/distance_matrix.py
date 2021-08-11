import numpy as np
from scipy.spatial.distance import cdist


def create_distance_matrix(arr: np.ndarray) -> np.ndarray:
    """
    This calculates a distance matrix D of nxn size for given array of n elements.
    The (i,j) element of matrix will have the euclidean distance from ith point to the jth point.
    """

    return cdist(arr, arr, metric="euclidean")


def find_m_nearest_idx(m: int, p_i: int, distance_matrix: np.ndarray, selected_idx: np.ndarray) -> np.ndarray:
    """
    Finds the indices of m nearest neighbors of point p given by p_i.
    (Notice that p_i is the index of p, not the actual point itself and this returns indices of the m nearest neighbors,
    not the points themselves.) Uses distance matrix D to calculate the nearest points. selected_idx can be used to
    filter only the required points. Only the indices provided in selected_idx array will be considered for finding
    neighbors.

    Note: p will not be considered as a nearest neighbor to p, even if p_i is provided in selected_idx.

    Warning: If the selected_idx filter is given less number of indices than m,
    this will throw an error return the selected_idx itself.

    :param m: Number of nearest neighbors to find.
    :param p_i: Index of point to find the neighbors.
    :param distance_matrix: Distance matrix.
    :param selected_idx: Filter for the required points.
    :return: M nearest neighbor indices.
    """

    if len(selected_idx) == 0:
        print("Warning: find_m_nearest_idx received 0 points to find neighbors.")
        return np.zeros(0)
    if len(selected_idx) < m:
        return selected_idx[:]

    dist_i = np.ones(len(distance_matrix)) * np.inf
    dist_i[selected_idx] = distance_matrix[p_i, selected_idx]
    dist_i[p_i] = np.inf  # this does not mutate original D
    sorted_dist_idx = np.argpartition(dist_i, m)
    return sorted_dist_idx[:m]

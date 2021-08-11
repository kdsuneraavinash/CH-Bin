import numpy as np
from scipy.spatial.distance import cdist

from bin_x.core.convex_hull.solve_qp import solve_qp


def distance_convex(query: np.ndarray, points: np.ndarray, solver: str = "quadprog") -> (float, np.ndarray, np.ndarray):
    """
    Finds distance to the convex hull using the given solver.

    :param query: Point to find the distance from.
    :param points: Points forming the target convex hull.
    :param solver: Solver to use.
    :return: A tuple containing the convex hull distance, the projection point coordinates
                and the optimized alpha value.
    """

    n = len(points)

    mat_a = np.ones((1, n))  # first equality constraint
    b = np.ones(1)  # first constraint RHS having one

    mat_g = -np.eye(n)  # second greater than zero constraint
    h = np.zeros(n)  # second constraint RHS having 0

    mat_x = points.copy().T  # points in cluster
    x = query.copy().T  # query point

    mat_p = 2 * np.matmul(mat_x.T, mat_x)
    q = (-2 * np.matmul(x.T, mat_x)).T

    alpha = solve_qp(mat_p, q, solver, mat_g, h, mat_a, b)
    proj = np.matmul(alpha, mat_x.T)
    return np.linalg.norm(proj - x), proj, alpha


def find_distance_matrix(arr: np.ndarray) -> np.ndarray:
    """
    This calculates a distance matrix D of nxn size for given array of n elements.
    The (i,j) element of matrix will have the euclidean distance from ith point to the jth point.
    """

    return cdist(arr, arr, metric="euclidean")


def find_m_nearest_idx(m: int, p_i: int, mat_dist: np.ndarray, selected_idx: np.ndarray) -> np.ndarray:
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
    :param mat_dist: Distance matrix.
    :param selected_idx: Filter for the required points.
    :return: M nearest neighbor indices.
    """

    if len(selected_idx) == 0:
        print("Warning: find_m_nearest_idx received 0 points to find neighbors.")
        return np.zeros(0)
    if len(selected_idx) < m:
        return selected_idx[:]

    dist_i = np.ones(len(mat_dist)) * np.inf
    dist_i[selected_idx] = mat_dist[p_i, selected_idx]
    dist_i[p_i] = np.inf  # this does not mutate original D
    sorted_dist_idx = np.argpartition(dist_i, m)
    return sorted_dist_idx[:m]

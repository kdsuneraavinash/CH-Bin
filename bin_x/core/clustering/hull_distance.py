import numpy as np

from bin_x.core.clustering.solve_qp import solve_qp


def convex_hull_distance(query: np.ndarray, points: np.ndarray, solver: str = "quadprog") -> float:
    """
    Finds distance to the convex hull using the given solver.

    :param query: Point to find the distance from.
    :param points: Points forming the target convex hull.
    :param solver: Solver to use. (quadprog/cvxopt)
    :return: The convex hull distance from the query point to the convex hull formed by the points given.
    """

    n = len(points)

    mat_a = np.ones((1, n))  # first equality constraint
    vec_b = np.ones(1)  # first constraint RHS having one

    mat_g = -np.eye(n)  # second greater than zero constraint
    vec_h = np.zeros(n)  # second constraint RHS having 0

    mat_x = points.copy().T  # points in cluster
    x = query.copy().T  # query point

    mat_p = 2 * np.matmul(mat_x.T, mat_x)
    vec_q = (-2 * np.matmul(x.T, mat_x)).T

    alpha = solve_qp(mat_p, vec_q, mat_g, vec_h, mat_a, vec_b, solver)
    proj = np.matmul(alpha, mat_x.T)
    return np.linalg.norm(proj - x)


def affine_hull_distance(query: np.ndarray, points: np.ndarray, solver: str = "quadprog") -> float:
    """
    Finds distance to the affine hull using the given solver.

    :param query: Point to find the distance from.
    :param points: Points forming the target affine hull.
    :param solver: Solver to use. (quadprog/cvxopt)
    :return: The affine hull distance from the query point to the affine hull formed by the points given.
    """

    raise NotImplementedError()

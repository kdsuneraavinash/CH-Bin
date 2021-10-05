import numpy as np
import scipy.linalg

from ch_bin.core.clustering.solve_qp import solve_qp


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

    mat_x = points.T  # points in cluster
    mat_x_t = points
    x = query.T  # query point
    x_t = query

    mat_p = 2 * np.matmul(mat_x_t, mat_x)
    vec_q = (-2 * np.matmul(x_t, mat_x)).T

    alpha = solve_qp(mat_p, vec_q, mat_g, vec_h, mat_a, vec_b, solver)
    proj = np.matmul(alpha, mat_x_t)
    return np.linalg.norm(proj - x)


def affine_hull_distance_qp(query: np.ndarray, points: np.ndarray, solver: str = "quadprog") -> float:
    """
    Finds distance to the affine hull using the given solver.

    :param query: Point to find the distance from.
    :param points: Points forming the target convex hull.
    :param solver: Solver to use. (quadprog/cvxopt)
    :return: The convex hull distance from the query point to the convex hull formed by the points given.
    """

    n = len(points)

    mat_a = np.ones((1, n))  # first equality constraint
    vec_b = np.ones(1)  # first constraint RHS having one

    mat_g = np.zeros(shape=(0, n))  # no second constraint
    vec_h = np.zeros(0)  # second constraint RHS having 0

    mat_x = points.T  # points in cluster
    mat_x_t = points
    x = query.T  # query point
    x_t = query

    mat_p = 2 * np.matmul(mat_x_t, mat_x)
    vec_q = (-2 * np.matmul(x_t, mat_x)).T

    alpha = solve_qp(mat_p, vec_q, mat_g, vec_h, mat_a, vec_b, solver)
    proj = np.matmul(alpha, mat_x_t)
    return np.linalg.norm(proj - x)


def affine_hull_distance(query: np.ndarray, points: np.ndarray) -> float:
    """
    Finds distance to the affine hull using the given solver.

    :param query: Point to find the distance from.
    :param points: Points forming the target affine hull.
    :return: The affine hull distance from the query point to the affine hull formed by the points given.
    """

    # Calculate the average each dimension
    mean_vec = points.mean(axis=0)
    relative_vecs = points - mean_vec
    # Find orthogonal basis
    basis = scipy.linalg.orth(relative_vecs.T)
    # Find the min distance
    ortho_proj_op = basis @ np.linalg.inv(basis.T @ basis) @ basis.T
    eye = np.eye(ortho_proj_op.shape[0])
    dist_vec = (eye - ortho_proj_op) @ (query - mean_vec)
    return np.linalg.norm(dist_vec)


def calculate_distance(x: np.ndarray, mat_p: np.ndarray, qp_solver: str, metric: str) -> float:
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
        return affine_hull_distance(x, mat_p)
    if metric == "affine-qp":
        return affine_hull_distance_qp(x, mat_p, solver=qp_solver)
    raise NotImplementedError(f"Metric {metric} not implemented")

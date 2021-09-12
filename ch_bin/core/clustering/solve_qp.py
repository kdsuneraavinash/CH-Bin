"""
Following functions will generalize quadprog and cvxopt tools to be applied using the same interface.

These functions will minimize (1/2) x^T P x + q^T x under G x <= h and A x = b. Here, x is the vector of optimization
variables x1 , ..., xn. The matrix P and vector q are used to define any quadratic objective function on these
variables, while the matrix-vector couples (G,h) and (A,b) respectively define inequality and equality constraints.

Derived from: https://scaron.info/blog/quadratic-programming-in-python.html
"""

import cvxopt
import numpy as np
import quadprog

from ch_bin.core.clustering.positive_def import nearest_positive_definite


def _quadprog_solve_qp(
    mat_p: np.ndarray,
    vec_q: np.ndarray,
    mat_g: np.ndarray,
    vec_h: np.ndarray,
    mat_a: np.ndarray,
    vec_b: np.ndarray,
) -> np.ndarray:
    """
    Use CVXOPT to calculate to solve the quadratic equation.

    quadprog requires a symmetric positive definite matrix.
    So, first, a closest positive-definite matrix is found.
    Generally, this matrix is very close to the original matrix.
    Thus, the result is good enough for the use case.

    This function will minimize (1/2) x^T P x + q^T x under G x <= h and A x = b.

    :param mat_p: Matrix P.
    :param vec_q: vector q.
    :param mat_g: Matrix G.
    :param vec_h: vector h.
    :param mat_a: Matrix A.
    :param vec_b: Vector b.
    :return: The solution vector.
    """
    qp_g = nearest_positive_definite(mat_p)

    qp_a = -vec_q
    qp_c = -np.concatenate((mat_a, mat_g), axis=0).T
    qp_b = -np.concatenate((vec_b, vec_h), axis=0)
    meq = mat_a.shape[0]

    return quadprog.solve_qp(qp_g, qp_a, qp_c, qp_b, meq)[0]


def _cvxopt_solve_qp(
    mat_p: np.ndarray,
    vec_q: np.ndarray,
    mat_g: np.ndarray,
    vec_h: np.ndarray,
    mat_a: np.ndarray,
    vec_b: np.ndarray,
) -> np.ndarray:
    """
    Use CVXOPT to calculate to solve the quadratic equation.

    CVXOPT assumes that you provide a symmetric cost matrix right away.
    CVXOPT won't check this, and will return a wrong results if it isn't.
    So the matrix is first transformed to a symmetric matrix.

    This function will minimize (1/2) x^T P x + q^T x under G x <= h and A x = b.

    :param mat_p: Matrix P.
    :param vec_q: vector q.
    :param mat_g: Matrix G.
    :param vec_h: vector h.
    :param mat_a: Matrix A.
    :param vec_b: Vector b.
    :return: The solution vector.
    """
    mat_p = 0.5 * (mat_p + mat_p.T)

    args = [
        cvxopt.matrix(mat_p),
        cvxopt.matrix(vec_q),
        cvxopt.matrix(mat_g),
        cvxopt.matrix(vec_h),
        cvxopt.matrix(mat_a),
        cvxopt.matrix(vec_b),
    ]
    cvxopt.solvers.options["show_progress"] = False
    sol = cvxopt.solvers.qp(*args)
    if "optimal" not in sol["status"]:
        raise ValueError("CVXOPT failed because sub-optimal solution.")
    return np.array(sol["x"]).reshape((mat_p.shape[1],))


def solve_qp(
    mat_p: np.ndarray,
    vec_q: np.ndarray,
    mat_g: np.ndarray,
    vec_h: np.ndarray,
    mat_a: np.ndarray,
    vec_b: np.ndarray,
    solver: str = "quadprog",
) -> np.ndarray:
    """
    Following function can be used to solve a quadratic equation using quadprog/cvxopt.

    Due to some reason, even if we convert the matrix G of quadprog to a positive definite matrix,
    the quadprog algorithm still complains with matrix G is not positive definite value error.
    cvxopt is not susceptible to this problem even if that algorithm is much slower than quadprog.
    So, as a trade-off, when the algorithm is set as quadprog, we still fallback to cvxopt if quadprog
    fails for some reason.

    This function will minimize (1/2) x^T P x + q^T x under G x <= h and A x = b.

    :param mat_p: Matrix P.
    :param vec_q: vector q.
    :param mat_g: Matrix G.
    :param vec_h: vector h.
    :param mat_a: Matrix A.
    :param vec_b: Vector b.
    :param solver: Solver algorithm to use. (quadprog/cvxopt)
    :return: The solution vector.
    """
    if solver == "quadprog":
        try:
            return _quadprog_solve_qp(mat_p, vec_q, mat_g, vec_h, mat_a, vec_b)
        except ValueError:
            return _cvxopt_solve_qp(mat_p, vec_q, mat_g, vec_h, mat_a, vec_b)
    elif solver == "cvxopt":
        return _cvxopt_solve_qp(mat_p, vec_q, mat_g, vec_h, mat_a, vec_b)
    raise NotImplementedError(f"Unknown solver {solver}")

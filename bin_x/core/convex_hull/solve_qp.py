"""
Following functions will generalize quadprog and cvxopt tools to be applied using the same interface.

These functions will minimize (1/2) x^T P x + q^T x under G x <= h and A x = b. Here, x is the vector of optimization
variables x1 , ..., xn. The matrix P and vector q are used to define any quadratic objective function on these
variables, while the matrix-vector couples (G,h) and (A,b) respectively define inequality and equality constraints.

Derived from: https://scaron.info/blog/quadratic-programming-in-python.html
"""
from typing import Optional

import cvxopt
import numpy as np
import quadprog

from bin_x.core.convex_hull.pos_def import nearest_positive_definite
from bin_x.core.utils import warn_once


def _quadprog_solve_qp(
    mat_p: np.ndarray,
    q: np.ndarray,
    mat_g: Optional[np.ndarray] = None,
    h: Optional[np.ndarray] = None,
    mat_a: Optional[np.ndarray] = None,
    b: Optional[np.ndarray] = None,
) -> np.array:
    # quadprog requires a symmetric positive definite matrix.
    # So we need to find closest pd matrix to G.
    qp_g = 0.5 * (mat_p + mat_p.T)
    qp_g = nearest_positive_definite(qp_g)

    qp_a = -q
    if mat_a is not None:
        qp_c = -np.vstack([mat_a, mat_g]).T
        qp_b = -np.hstack([b, h])
        meq = mat_a.shape[0]
    else:
        qp_c = -mat_g.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_g, qp_a, qp_c, qp_b, meq)[0]


def _cvxopt_solve_qp(
    mat_p: np.ndarray,
    q: np.ndarray,
    mat_g: Optional[np.ndarray] = None,
    h: Optional[np.ndarray] = None,
    mat_a: Optional[np.ndarray] = None,
    b: Optional[np.ndarray] = None,
) -> np.ndarray:
    # CVXOPT assumes that you provide a symmetric cost matrix right away.
    # CVXOPT won't check this, and will return a wrong results if it isnt.
    mat_p = 0.5 * (mat_p + mat_p.T)

    args = [cvxopt.matrix(mat_p), cvxopt.matrix(q)]
    if mat_g is not None:
        args.extend([cvxopt.matrix(mat_g), cvxopt.matrix(h)])
        if mat_a is not None:
            args.extend([cvxopt.matrix(mat_a), cvxopt.matrix(b)])
    cvxopt.solvers.options["show_progress"] = False
    sol = cvxopt.solvers.qp(*args)
    if "optimal" not in sol["status"]:
        raise ValueError("CVXOPT failed because sub-optimal solution.")
    return np.array(sol["x"]).reshape((mat_p.shape[1],))


def solve_qp(
    mat_p: np.ndarray,
    q: np.ndarray,
    solver: str,
    mat_g: Optional[np.ndarray] = None,
    h: Optional[np.ndarray] = None,
    mat_a: Optional[np.ndarray] = None,
    b: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Following function can be used to solve a quadratic equation using quadprog/cvxopt.

    Due to some reason, even if we convert the matrix G of quadprog to a positive definite matrix,
    the quadprog algorithm still complains with matrix G is not positive definite value error.
    cvxopt is not susceptible to this problem even if that algorithm is much slower than quadprog.
    So, as a trade-off, when the algorithm is set as quadprog, we still fallback to cvxopt if quadprog
    fails for some reason.
    """
    if solver == "quadprog":
        try:
            return _quadprog_solve_qp(mat_p, q, mat_g, h, mat_a, b)
        except ValueError as e:
            warn_once("FBK_QUADPROG", f"Falling back to cvxopt. ({e})")
            return _cvxopt_solve_qp(mat_p, q, mat_g, h, mat_a, b)
    elif solver == "cvxopt":
        return _cvxopt_solve_qp(mat_p, q, mat_g, h, mat_a, b)
    raise NotImplementedError(f"Unknown solver {solver}")

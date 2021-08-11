"""
Module containing functions relating to positive definite matrices.

Quadratic programming using quadprog requires calculating a positive definite matrix for a given matrix A.
This solution provides the implementation to find the nearest positive definite matrix.
Implementation derived from
https://stackoverflow.com/questions/43238173/python-convert-matrix-to-positive-semi-definite/43244194#43244194

This is a Python/Numpy port of John D'Errico's nearestSPD MATLAB code [1], which credits [2].
However, note that spacing calculation is different from [1].
MATLAB's chol Cholesky decomposition will accept matrices with exactly 0-eigenvalue, whereas Numpy's will not.
So where [1] uses eps(mineig) (where eps is Matlab for np.spacing), we use the above definition.
CAVEAT: our spacing will be much larger than [1]'s eps(mineig), since mineig is usually on the order of 1e-16,
and eps(1e-16) is on the order of 1e-34, whereas spacing will, for Gaussian random matrices of small dimension,
be on the order of 1e-16. In practice, both ways converge.

[1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd \n
[2] N.J. Higham, "Computing a nearest symmetric positive semidefinite matrix" (1988):
https://doi.org/10.1016/0024-3795(88)90223-6
"""

import numpy as np


def nearest_positive_definite(mat_a: np.ndarray) -> np.ndarray:
    """
    Calculates the nearest positive definite matrix to the given matrix.

    :param mat_a: Input matrix.
    :return: Same matrix but with small adjustments to make it positive definite.
    """

    mat_b = (mat_a + mat_a.T) / 2
    _, s, mat_v = np.linalg.svd(mat_b)
    mat_h = np.dot(mat_v.T, np.dot(np.diag(s), mat_v))
    mat_a2 = (mat_b + mat_h) / 2
    mat_a3 = (mat_a2 + mat_a2.T) / 2
    if is_positive_definite(mat_a3):
        return mat_a3
    spacing = np.spacing(np.linalg.norm(mat_a))
    mat_i = np.eye(mat_a.shape[0])
    k = 1
    while not is_positive_definite(mat_a3):
        mineig = np.min(np.real(np.linalg.eigvals(mat_a3)))
        mat_a3 += mat_i * (-mineig * k ** 2 + spacing)
        k += 1
    return mat_a3


def is_positive_definite(mat_a) -> bool:
    """
    Returns whether the given matrix is positive definite.

    :param mat_a: Input matrix.
    :return: Whether the matrix is positive definite.
    """

    try:
        _ = np.linalg.cholesky(mat_a)
        return True
    except np.linalg.LinAlgError:
        return False

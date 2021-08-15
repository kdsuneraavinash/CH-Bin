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

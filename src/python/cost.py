from typing import Union

import numpy as np
import scipy.sparse


def calcCost(
    pos: np.ndarray, A: Union[np.ndarray, scipy.sparse.csr_matrix], k: float
) -> float:
    shape = A.shape
    assert shape is not None
    nrows, ncols = shape
    assert nrows is not None and ncols is not None
    assert pos.shape[0] == nrows == ncols
    delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
    distance = np.linalg.norm(delta, axis=-1)
    EPS = 1e-10
    if isinstance(A, np.ndarray):
        return np.sum(A * distance**3 / (3 * k) - (k**2) * np.log(distance + EPS))
    if scipy.sparse.issparse(A):
        return np.sum(
            A.multiply(distance**3) / (3 * k) - (k**2) * np.log(distance + EPS)
        )
    raise TypeError("A must be either np.ndarray or scipy.sparse.csr_matrix")


if __name__ == "__main__":
    pos = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    A = np.array([[1.0, 0.5, 0.5], [0.5, 1.0, 0.5], [0.5, 0.5, 1.0]])
    print(calcCost(pos, A, 1))
    A = scipy.sparse.csr_matrix(A)
    print(calcCost(pos, A, 1))

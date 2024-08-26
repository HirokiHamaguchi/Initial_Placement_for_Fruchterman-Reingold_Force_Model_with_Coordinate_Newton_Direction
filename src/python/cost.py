import numpy as np
import scipy.sparse
from typing import Union


def calcCost(
    pos: np.ndarray, A: Union[np.ndarray, scipy.sparse.csr_matrix], k: float
) -> float:
    assert pos.shape[0] == A.shape[0] == A.shape[1]
    delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
    distance = np.linalg.norm(delta, axis=-1)
    EPS = 1e-10
    if type(A) == np.ndarray:
        return np.sum(A * distance**3 / (3 * k) - (k**2) * np.log(distance + EPS))
    else:
        return np.sum(
            A.multiply(distance**3) / (3 * k) - (k**2) * np.log(distance + EPS)
        )


if __name__ == "__main__":
    pos = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    A = np.array([[1.0, 0.5, 0.5], [0.5, 1.0, 0.5], [0.5, 0.5, 1.0]])
    print(calcCost(pos, A, 1))
    A = scipy.sparse.csr_matrix(A)
    print(calcCost(pos, A, 1))

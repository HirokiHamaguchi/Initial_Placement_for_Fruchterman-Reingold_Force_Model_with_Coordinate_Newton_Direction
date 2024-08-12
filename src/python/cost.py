import numpy as np


def cost(pos: np.ndarray, A: np.ndarray, k: float) -> float:
    assert pos.shape[0] == A.shape[0] == A.shape[1]
    delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
    distance = np.linalg.norm(delta, axis=-1)
    EPS = 1e-10
    return np.sum(A * (distance**3) / (3 * k) - (k**2) * np.log(distance + EPS))


if __name__ == "__main__":
    pos = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    A = np.array([[1.0, 0.5, 0.5], [0.5, 1.0, 0.5], [0.5, 0.5, 1.0]])
    print(cost(pos, A, 1))

from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse


def visMat(mat: Union[np.ndarray, scipy.sparse.coo_matrix]) -> None:
    plt.figure(figsize=(10, 10))
    if scipy.sparse.issparse(mat):
        mat_sparse = scipy.sparse.coo_matrix(mat)
        matVis = mat_sparse.toarray().astype(float)
    else:
        matVis = mat.copy().astype(float)
    matVis[matVis == 0] = np.nan
    plt.imshow(matVis, cmap="viridis", interpolation="none")
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    visMat(np.array([[1, 2, 3], [2, 1, 4], [3, 4, 1]]))
    visMat(scipy.sparse.coo_matrix([[1, 2, 3], [2, 1, 4], [3, 4, 1]]))

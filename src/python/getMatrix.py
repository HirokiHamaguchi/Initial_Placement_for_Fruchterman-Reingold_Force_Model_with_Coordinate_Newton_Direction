import glob
import os
from typing import Generator, Union

import numpy as np
import scipy.io
import scipy.sparse


def getMatrixByName(
    name: str, asDense: bool = False
) -> Union[np.ndarray, scipy.sparse.coo_matrix]:
    assert os.getcwd().count("FruchtermanReingoldByRandomSubspace") == 1
    path = (
        os.getcwd().split("FruchtermanReingoldByRandomSubspace")[0]
        + "FruchtermanReingoldByRandomSubspace/data/"
    )
    path = path + name + ".mtx"
    assert os.path.exists(path), f"Matrix({name}) not found"
    matrix = scipy.io.mmread(path)
    if scipy.sparse.issparse(matrix):
        matrix_sparse = scipy.sparse.coo_matrix(matrix)
        if asDense:
            return matrix_sparse.toarray()
        return matrix_sparse
    matrix_dense = np.asarray(matrix)
    return matrix_dense


def getMatrixes(
    asDense: bool = False,
) -> Generator[Union[np.ndarray, scipy.sparse.coo_matrix], None, None]:
    assert os.path.exists("../../data"), "Data folder not found"
    names = glob.glob("../../data/*.mtx")
    for name in names:
        yield getMatrixByName(name[11:-4], asDense)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print(getMatrixByName("jagmesh1"))
    for mat in getMatrixes():
        shape = mat.shape
        assert shape is not None
        print(shape)
        assert shape[0] == shape[1]

import os
import glob
import networkx as nx
import numpy as np
import scipy.sparse

from src.python.getMatrix import getMatrixByName
from src.python.vis.visGraph import visGraph
from src.python.spring_layout import spring_layout
from packaging import version
import PIL


def visAnimate(matrixName: str, methodName: str):
    mat = getMatrixByName(matrixName)
    if scipy.sparse.issparse(mat):
        mat.setdiag(0)
        mat.eliminate_zeros()
        mat.data = np.abs(mat.data)
    else:
        mat[np.diag_indices_from(mat)] = 0
        mat.data = np.abs(mat.data)

    assert os.path.exists("temp"), "temp directory does not exist"

    G = nx.Graph(mat)
    for i, pos in enumerate(spring_layout(G, method=methodName, iterations=150)):
        print(f"iteration: {i}", end="\r")
        visGraph(
            G,
            pos,
            title=f"{matrixName}_{methodName}_{i}",
            savePath=f"temp/{i:03}.png",
        )
    print()

    if version.parse(PIL.__version__) < version.parse("3.4"):
        # https://stackoverflow.com/questions/24688802/saving-an-animated-gif-in-pillow
        raise Exception(
            "Pillow in version not supporting making animated GIFs",
            "you need to upgrade library version",
            "see release notes in",
            "https://pillow.readthedocs.io/en/latest/releasenotes/3.4.0.html#append-images-to-gif",
        )
    else:
        files = sorted(glob.glob(f"temp/*.png"))
        images = [PIL.Image.open(file) for file in files]
        images[0].save(
            f"{matrixName}_{methodName}.gif",
            save_all=True,
            append_images=images[1:],
            duration=100,
            loop=0,
        )


if __name__ == "__main__":
    visAnimate("jagmesh1", methodName="BFGS")

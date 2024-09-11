import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from typing import Union


def visGraph(
    G: Union[nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph],
    pos: dict,
    _dirs: np.ndarray = None,
    title: str = None,
    savePath: str = None,
    node_size: int = 50,
    width: float = 1.0,
) -> None:
    cmap = plt.get_cmap("jet")
    n = G.number_of_nodes()
    colorMap = np.array([cmap(i / (n - 1)) for i in range(n)])
    nx.draw(G, pos, node_size=node_size, node_color=colorMap, width=width)
    if _dirs is not None and np.any(_dirs != 0.0):
        dirs = _dirs.copy()
        assert pos.shape == (n, 2)
        assert dirs.shape == (n, 2)
        dirs /= np.max(np.linalg.norm(dirs, axis=1)) / 0.3
        for i in range(n):
            plt.arrow(
                pos[i][0],
                pos[i][1],
                dirs[i, 0],
                dirs[i, 1],
                head_width=0.05,
                head_length=0.05,
                fc="k",
                ec="k",
            )
    if title:
        plt.title(title)
    if savePath:
        plt.savefig(savePath, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    mat = np.array([[1, 2, 3], [2, 1, 4], [3, 4, 1]])
    G = nx.Graph(mat)
    pos = nx.spring_layout(G)
    visGraph(G, pos)
